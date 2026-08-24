from __future__ import annotations
import os
import numpy as np
import pandas as pd
import biom
import skbio

try:
    import qiime2
except ImportError:
    class DummyQIIME2:
        Metadata = object
    qiime2 = DummyQIIME2()

from scipy.spatial.distance import pdist, squareform
from concurrent.futures import ProcessPoolExecutor

def _run_bnti_perm(args):
    seed, n_asvs, n_samples, P_val, D_matrix = args
    np.random.seed(seed)
    perm = np.random.permutation(n_asvs)
    D_perm = D_matrix[perm, :][:, perm]
    M_min = np.zeros((n_asvs, n_samples))
    for j in range(n_samples):
        present_j = P_val[j, :] > 0
        if np.any(present_j):
            M_min[:, j] = np.min(D_perm[:, present_j], axis=1)
    mntd = P_val @ M_min
    return 0.5 * (mntd + mntd.T)

def _run_rc_perm(args):
    seed, n_samples, seq_depths, regional_prob = args
    np.random.seed(seed)
    X_null = np.array([
        np.random.multinomial(seq_depths[i], regional_prob)
        for i in range(n_samples)
    ], dtype=float)
    sums = X_null.sum(axis=1, keepdims=True)
    sums[sums == 0] = 1.0
    rel_null = X_null / sums
    return squareform(pdist(rel_null, metric='braycurtis'))

def calculate_processes(
    table: biom.Table,
    tree: skbio.TreeNode = None,
    metadata: qiime2.Metadata = None,
    column: str = None,
    permutations: int = 999,
    min_frequency: int = 0,
    n_jobs: int = 1
) -> biom.Table:
    """
    QIIME 2 Method Action: Quantifies the 5 Stegen ecological assembly processes
    (Homogeneous Selection, Variable Selection, Dispersal Limitation, Homogenizing Dispersal, Undominated Drift)
    and Ning's Normalized Stochasticity Ratio (NST) using phylogenetic null models across experimental groups.
    """
    if table.is_empty() or table.matrix_data.sum() == 0:
        raise ValueError("Provided FeatureTable[Frequency] is empty.")

    if tree is None:
        raise ValueError("A rooted phylogenetic tree (Phylogeny[Rooted]) is required for BetaNTI null model calculation.")

    table_df = pd.DataFrame(
        table.matrix_data.toarray().T,
        index=table.ids(axis='sample'),
        columns=table.ids(axis='observation')
    )

    sample_depths = table_df.sum(axis=1)
    if len(sample_depths) > 0 and sample_depths.min() > 0:
        depth_ratio = sample_depths.max() / float(sample_depths.min())
        if depth_ratio > 3.0:
            import warnings
            warnings.warn(
                f"Sequencing depth varies by {depth_ratio:.1f}x across samples. "
                "It is strongly recommended to rarefy your feature table (`qiime feature-table rarefy`) "
                "prior to running q2-eco-processes to avoid library-size bias in RCbray null models.",
                UserWarning
            )

    if min_frequency > 0:
        sums = table_df.sum(axis=0)
        table_df = table_df.loc[:, sums >= min_frequency]

    tree_tips = set(t.name for t in tree.tips())
    common_asvs = [a for a in table_df.columns if a in tree_tips]

    if len(common_asvs) < 2:
        raise ValueError("Fewer than 2 matching ASVs found between FeatureTable and Phylogeny[Rooted].")

    table_df = table_df[common_asvs]
    sample_ids = list(table_df.index)
    n_samples = len(sample_ids)
    n_asvs = len(common_asvs)

    if n_samples < 2:
        raise ValueError("At least 2 samples are required to calculate ecological assembly processes.")

    dm = tree.tip_tip_distances()
    dm_df = pd.DataFrame(dm.data, index=dm.ids, columns=dm.ids)
    D_matrix = dm_df.loc[common_asvs, common_asvs].values.astype(float)
    P_val = table_df.div(table_df.sum(axis=1), axis=0).fillna(0.0).values.astype(float)

    # 1. Observed BetaMNTD
    M_min_obs = np.zeros((n_asvs, n_samples))
    for j in range(n_samples):
        present_j = P_val[j, :] > 0
        if np.any(present_j):
            M_min_obs[:, j] = np.min(D_matrix[:, present_j], axis=1)

    mntd_obs_dir = P_val @ M_min_obs
    beta_mntd_obs = 0.5 * (mntd_obs_dir + mntd_obs_dir.T)

    workers = min(n_jobs if n_jobs > 0 else os.cpu_count(), os.cpu_count() or 1)
    bnti_args = [
        (s, n_asvs, n_samples, P_val, D_matrix)
        for s in range(42, 42 + permutations)
    ]

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            null_mntds = list(executor.map(_run_bnti_perm, bnti_args))
    else:
        null_mntds = [_run_bnti_perm(arg) for arg in bnti_args]

    null_mntds_arr = np.array(null_mntds)
    mean_null_mntd = np.mean(null_mntds_arr, axis=0)
    std_null_mntd = np.std(null_mntds_arr, axis=0)
    std_null_mntd[std_null_mntd == 0] = 1.0

    beta_nti = (beta_mntd_obs - mean_null_mntd) / std_null_mntd
    np.fill_diagonal(beta_nti, 0.0)

    # 2. Observed Bray-Curtis & Raup-Crick (RCbray) null model
    bc_obs = squareform(pdist(P_val, metric='braycurtis'))
    seq_depths = table_df.sum(axis=1).values.astype(int)
    seq_depths[seq_depths <= 0] = 1

    regional_prob = table_df.sum(axis=0).values.astype(float)
    reg_sum = regional_prob.sum()
    if reg_sum > 0:
        regional_prob /= reg_sum
    else:
        regional_prob = np.ones(n_asvs) / n_asvs

    rc_args = [
        (s, n_samples, seq_depths, regional_prob)
        for s in range(1000, 1000 + permutations)
    ]

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            null_bcs = list(executor.map(_run_rc_perm, rc_args))
    else:
        null_bcs = [_run_rc_perm(arg) for arg in rc_args]

    null_bcs_arr = np.array(null_bcs)
    rc_bray = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            if i == j:
                continue
            obs_val = bc_obs[i, j]
            null_vals = null_bcs_arr[:, i, j]
            count_gt = np.sum(null_vals > obs_val)
            count_eq = np.sum(null_vals == obs_val)
            rc_bray[i, j] = 2.0 * ((count_gt + 0.5 * count_eq) / float(permutations)) - 1.0

    processes = [
        "Homogeneous_Selection",
        "Variable_Selection",
        "Dispersal_Limitation",
        "Homogenizing_Dispersal",
        "Undominated_Drift",
        "Normalized_Stochasticity_Ratio",
        "Mean_BetaNTI",
        "Mean_RCbray"
    ]

    if metadata is not None:
        meta_df = metadata.to_dataframe() if hasattr(metadata, 'to_dataframe') else metadata
    else:
        meta_df = None

    if meta_df is not None and column is not None and column in meta_df.columns:
        valid_meta = meta_df.loc[meta_df.index.intersection(sample_ids)]
        groups = list(valid_meta[column].unique())
        group_counts = pd.DataFrame(0.0, index=groups, columns=processes)

        for grp in groups:
            grp_samples = valid_meta[valid_meta[column] == grp].index
            indices = [sample_ids.index(s) for s in grp_samples if s in sample_ids]

            bnti_list = []
            rc_list = []

            for i in indices:
                for j in indices:
                    if i == j:
                        continue
                    b_nti = beta_nti[i, j]
                    rc = rc_bray[i, j]
                    bnti_list.append(b_nti)
                    rc_list.append(rc)

                    if b_nti < -2.0:
                        group_counts.loc[grp, "Homogeneous_Selection"] += 1.0
                    elif b_nti > 2.0:
                        group_counts.loc[grp, "Variable_Selection"] += 1.0
                    elif rc > 0.95:
                        group_counts.loc[grp, "Dispersal_Limitation"] += 1.0
                    elif rc < -0.95:
                        group_counts.loc[grp, "Homogenizing_Dispersal"] += 1.0
                    else:
                        group_counts.loc[grp, "Undominated_Drift"] += 1.0

            group_counts.loc[grp, "Mean_BetaNTI"] = np.mean(bnti_list) if len(bnti_list) > 0 else 0.0
            group_counts.loc[grp, "Mean_RCbray"] = np.mean(rc_list) if len(rc_list) > 0 else 0.0

        stoch_cols = ["Dispersal_Limitation", "Homogenizing_Dispersal", "Undominated_Drift"]
        det_cols = ["Homogeneous_Selection", "Variable_Selection"]
        total_counts = group_counts[stoch_cols + det_cols].sum(axis=1)
        stoch_counts = group_counts[stoch_cols].sum(axis=1)
        nst_val = (stoch_counts / total_counts.replace(0, 1.0) * 100.0).fillna(50.0)

        sums = total_counts.replace(0, 1.0)
        process_pct = group_counts[stoch_cols + det_cols].div(sums, axis=0) * 100.0
        process_pct["Normalized_Stochasticity_Ratio"] = nst_val
        process_pct["Mean_BetaNTI"] = group_counts["Mean_BetaNTI"]
        process_pct["Mean_RCbray"] = group_counts["Mean_RCbray"]
    else:
        sample_counts = pd.DataFrame(0.0, index=sample_ids, columns=processes)
        for i in range(n_samples):
            bnti_list = []
            rc_list = []
            for j in range(n_samples):
                if i == j:
                    continue
                b_nti = beta_nti[i, j]
                rc = rc_bray[i, j]
                bnti_list.append(b_nti)
                rc_list.append(rc)

                if b_nti < -2.0:
                    sample_counts.loc[sample_ids[i], "Homogeneous_Selection"] += 1.0
                elif b_nti > 2.0:
                    sample_counts.loc[sample_ids[i], "Variable_Selection"] += 1.0
                elif rc > 0.95:
                    sample_counts.loc[sample_ids[i], "Dispersal_Limitation"] += 1.0
                elif rc < -0.95:
                    sample_counts.loc[sample_ids[i], "Homogenizing_Dispersal"] += 1.0
                else:
                    sample_counts.loc[sample_ids[i], "Undominated_Drift"] += 1.0

            sample_counts.loc[sample_ids[i], "Mean_BetaNTI"] = np.mean(bnti_list) if len(bnti_list) > 0 else 0.0
            sample_counts.loc[sample_ids[i], "Mean_RCbray"] = np.mean(rc_list) if len(rc_list) > 0 else 0.0

        stoch_cols = ["Dispersal_Limitation", "Homogenizing_Dispersal", "Undominated_Drift"]
        det_cols = ["Homogeneous_Selection", "Variable_Selection"]
        total_counts = sample_counts[stoch_cols + det_cols].sum(axis=1)
        stoch_counts = sample_counts[stoch_cols].sum(axis=1)
        nst_val = (stoch_counts / total_counts.replace(0, 1.0) * 100.0).fillna(50.0)

        sums = total_counts.replace(0, 1.0)
        process_pct = sample_counts[stoch_cols + det_cols].div(sums, axis=0) * 100.0
        process_pct["Normalized_Stochasticity_Ratio"] = nst_val
        process_pct["Mean_BetaNTI"] = sample_counts["Mean_BetaNTI"]
        process_pct["Mean_RCbray"] = sample_counts["Mean_RCbray"]

    process_biom = biom.Table(
        process_pct.values.T,
        observation_ids=list(process_pct.columns),
        sample_ids=list(process_pct.index)
    )
    return process_biom

def calculate_matrices(
    table: biom.Table,
    tree: skbio.TreeNode = None,
    permutations: int = 999,
    min_frequency: int = 0,
    n_jobs: int = 1
) -> (skbio.DistanceMatrix, skbio.DistanceMatrix):
    """
    QIIME 2 Method Action: Calculates pairwise BetaNTI and Raup-Crick (RCbray) score matrices.
    """
    if table.is_empty() or table.matrix_data.sum() == 0:
        raise ValueError("Provided FeatureTable[Frequency] is empty.")

    if tree is None:
        raise ValueError("A rooted phylogenetic tree (Phylogeny[Rooted]) is required for BetaNTI calculation.")

    table_df = pd.DataFrame(
        table.matrix_data.toarray().T,
        index=table.ids(axis='sample'),
        columns=table.ids(axis='observation')
    )

    if min_frequency > 0:
        sums = table_df.sum(axis=0)
        table_df = table_df.loc[:, sums >= min_frequency]

    tree_tips = set(t.name for t in tree.tips())
    common_asvs = [a for a in table_df.columns if a in tree_tips]

    if len(common_asvs) < 2:
        raise ValueError("Fewer than 2 matching ASVs found between FeatureTable and Phylogeny[Rooted].")

    table_df = table_df[common_asvs]
    sample_ids = list(table_df.index)
    n_samples = len(sample_ids)
    n_asvs = len(common_asvs)

    if n_samples < 2:
        raise ValueError("At least 2 samples are required to calculate ecological assembly process matrices.")

    dm = tree.tip_tip_distances()
    dm_df = pd.DataFrame(dm.data, index=dm.ids, columns=dm.ids)
    D_matrix = dm_df.loc[common_asvs, common_asvs].values.astype(float)
    P_val = table_df.div(table_df.sum(axis=1), axis=0).fillna(0.0).values.astype(float)

    # 1. BetaNTI
    M_min_obs = np.zeros((n_asvs, n_samples))
    for j in range(n_samples):
        present_j = P_val[j, :] > 0
        if np.any(present_j):
            M_min_obs[:, j] = np.min(D_matrix[:, present_j], axis=1)

    mntd_obs_dir = P_val @ M_min_obs
    beta_mntd_obs = 0.5 * (mntd_obs_dir + mntd_obs_dir.T)

    workers = min(n_jobs if n_jobs > 0 else os.cpu_count(), os.cpu_count() or 1)
    bnti_args = [
        (s, n_asvs, n_samples, P_val, D_matrix)
        for s in range(42, 42 + permutations)
    ]

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            null_mntds = list(executor.map(_run_bnti_perm, bnti_args))
    else:
        null_mntds = [_run_bnti_perm(arg) for arg in bnti_args]

    null_mntds_arr = np.array(null_mntds)
    mean_null_mntd = np.mean(null_mntds_arr, axis=0)
    std_null_mntd = np.std(null_mntds_arr, axis=0)
    std_null_mntd[std_null_mntd == 0] = 1.0

    beta_nti = (beta_mntd_obs - mean_null_mntd) / std_null_mntd
    np.fill_diagonal(beta_nti, 0.0)

    # 2. RCbray
    bc_obs = squareform(pdist(P_val, metric='braycurtis'))
    seq_depths = table_df.sum(axis=1).values.astype(int)
    seq_depths[seq_depths <= 0] = 1

    regional_prob = table_df.sum(axis=0).values.astype(float)
    reg_sum = regional_prob.sum()
    if reg_sum > 0:
        regional_prob /= reg_sum
    else:
        regional_prob = np.ones(n_asvs) / n_asvs

    rc_args = [
        (s, n_samples, seq_depths, regional_prob)
        for s in range(1000, 1000 + permutations)
    ]

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            null_bcs = list(executor.map(_run_rc_perm, rc_args))
    else:
        null_bcs = [_run_rc_perm(arg) for arg in rc_args]

    null_bcs_arr = np.array(null_bcs)
    rc_bray = np.zeros((n_samples, n_samples))
    for i in range(n_samples):
        for j in range(n_samples):
            if i == j:
                continue
            obs_val = bc_obs[i, j]
            null_vals = null_bcs_arr[:, i, j]
            count_gt = np.sum(null_vals > obs_val)
            count_eq = np.sum(null_vals == obs_val)
            rc_bray[i, j] = 2.0 * ((count_gt + 0.5 * count_eq) / float(permutations)) - 1.0

    beta_nti_dm = skbio.DistanceMatrix(beta_nti, ids=sample_ids)
    rc_bray_dm = skbio.DistanceMatrix(rc_bray, ids=sample_ids)

    return beta_nti_dm, rc_bray_dm
