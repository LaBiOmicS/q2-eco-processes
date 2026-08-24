import os
import json
import numpy as np
import pandas as pd
import biom
import skbio
try:
    import qiime2
except ImportError:
    qiime2 = None
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import mantel, DistanceMatrix

def check_phylogenetic_signal(
    output_dir: str,
    table: biom.Table,
    tree: skbio.TreeNode,
    metadata: qiime2.Metadata,
    column: str,
    permutations: int = 999
) -> None:
    """
    QIIME 2 Visualizer Action: Evaluates the presence of phylogenetic signal by testing
    the correlation between environmental niche distance and phylogenetic distance using Mantel tests.
    """
    if table.is_empty():
        raise ValueError("Provided FeatureTable[Frequency] is empty.")

    if tree is None:
        raise ValueError("A rooted phylogenetic tree (Phylogeny[Rooted]) is required for phylogenetic signal analysis.")

    # Convert biom table to DataFrame
    table_df = pd.DataFrame(
        table.matrix_data.toarray().T,
        index=table.ids(axis='sample'),
        columns=table.ids(axis='observation')
    )

    # Filter metadata
    meta_df = metadata.to_dataframe() if hasattr(metadata, 'to_dataframe') else metadata
    if column not in meta_df.columns:
        raise ValueError(f"Column '{column}' not found in provided metadata.")

    # Ensure column is numeric
    try:
        env_series = pd.to_numeric(meta_df[column], errors='coerce').dropna()
    except Exception as e:
        raise ValueError(f"Column '{column}' must contain numeric environmental values: {e}")

    common_samples = list(table_df.index.intersection(env_series.index))
    if len(common_samples) < 3:
        raise ValueError(f"At least 3 valid samples with numeric '{column}' metadata are required for Mantel test.")

    table_df = table_df.loc[common_samples]
    env_series = env_series.loc[common_samples]

    # Align ASVs with tree tips
    tree_tips = set(t.name for t in tree.tips())
    common_asvs = [a for a in table_df.columns if a in tree_tips]
    if len(common_asvs) < 3:
        raise ValueError("At least 3 matching ASVs found between FeatureTable and Phylogeny[Rooted] are required.")

    table_df = table_df[common_asvs]

    # Calculate ASV environmental niche optima (abundance-weighted mean of environmental variable)
    rel_table = table_df.div(table_df.sum(axis=1), axis=0).fillna(0.0)
    asv_niche_optima = (rel_table.T @ env_series.values) / rel_table.sum(axis=0).values.clip(min=1e-9)

    # Niche distance matrix between ASVs
    niche_dist = squareform(pdist(asv_niche_optima.values.reshape(-1, 1), metric='euclidean'))
    niche_dm = DistanceMatrix(niche_dist, ids=common_asvs)

    # Phylogenetic distance matrix between ASVs
    phylo_dm_all = tree.tip_tip_distances()
    phylo_df = pd.DataFrame(phylo_dm_all.data, index=phylo_dm_all.ids, columns=phylo_dm_all.ids)
    phylo_dist = phylo_df.loc[common_asvs, common_asvs].values
    phylo_dm = DistanceMatrix(phylo_dist, ids=common_asvs)

    # Run Mantel Test
    r_stat, p_val, n_pairs = mantel(niche_dm, phylo_dm, permutations=permutations, method='pearson')
    r_stat_clean = round(float(r_stat), 4)
    p_val_clean = round(float(p_val), 4)

    # Prepare scatter plot data (Phylogenetic Distance vs Niche Distance)
    # Downsample points if ASV pairs > 2000 to keep HTML lightweight
    phylo_upper = phylo_dist[np.triu_indices_from(phylo_dist, k=1)]
    niche_upper = niche_dist[np.triu_indices_from(niche_dist, k=1)]

    if len(phylo_upper) > 2000:
        idx = np.random.choice(len(phylo_upper), size=2000, replace=False)
        x_pts = phylo_upper[idx].round(4).tolist()
        y_pts = niche_upper[idx].round(4).tolist()
    else:
        x_pts = phylo_upper.round(4).tolist()
        y_pts = niche_upper.round(4).tolist()

    # Determine signal status
    if p_val_clean < 0.05 and r_stat_clean > 0:
        signal_badge = "bg-success"
        signal_text = "SIGNIFICANT PHYLOGENETIC SIGNAL DETECTED"
        signal_desc = f"Phylogenetic distance correlates positively with environmental niche distance (Mantel r = {r_stat_clean}, p = {p_val_clean}). The assumption of Stegen's selection null models (BetaNTI) is scientifically validated for this dataset."
    else:
        signal_badge = "bg-warning text-dark"
        signal_text = "WEAK OR NON-SIGNIFICANT PHYLOGENETIC SIGNAL"
        signal_desc = f"No strong linear correlation was observed between phylogenetic distance and environmental niche distance (Mantel r = {r_stat_clean}, p = {p_val_clean}). Interpret BetaNTI Selection categories with caution."

    index_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>q2-eco-processes: Phylogenetic Signal Diagnostic</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <style>
        body {{ background-color: #f8f9fa; font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }}
        .header-box {{ background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); color: white; padding: 2rem; border-radius: 12px; margin-bottom: 2rem; }}
        .card-custom {{ border: none; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.05); margin-bottom: 2rem; background: white; }}
    </style>
</head>
<body>
<div class="container-fluid py-4 px-4">
    <div class="header-box">
        <h1 class="display-5 fw-bold mb-1">🧬 Phylogenetic Signal Analysis</h1>
        <p class="lead mb-0">Testing Ecological Niche Conservatism vs. Phylogenetic Distance (Mantel Test)</p>
    </div>

    <div class="row text-center mb-4">
        <div class="col-md-3">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">Environmental Variable</h6>
                <h3 class="fw-bold text-dark"><code>{column}</code></h3>
            </div>
        </div>
        <div class="col-md-3">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">Mantel Statistic (r)</h6>
                <h3 class="fw-bold text-primary">{r_stat_clean}</h3>
            </div>
        </div>
        <div class="col-md-3">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">p-value</h6>
                <h3 class="fw-bold text-info">{p_val_clean}</h3>
            </div>
        </div>
        <div class="col-md-3">
            <div class="card card-custom p-3">
                <h6 class="text-muted text-uppercase">Permutations</h6>
                <h3 class="fw-bold text-secondary">{permutations:,}</h3>
            </div>
        </div>
    </div>

    <div class="alert {signal_badge} p-4 mb-4 rounded-3 shadow-sm">
        <h4 class="fw-bold mb-2">📌 Status: {signal_text}</h4>
        <p class="m-0 fs-6">{signal_desc}</p>
    </div>

    <div class="card card-custom p-3 mb-4">
        <h5 class="fw-bold mb-3">📈 Phylogenetic Distance vs Environmental Niche Distance</h5>
        <div id="mantelScatter" style="height: 450px;"></div>
    </div>

    <div class="footer text-center text-muted py-3">
        <small>Generated with <strong>q2-eco-processes</strong> | LaBiOmics / UMC Microbiome Suite</small>
    </div>
</div>

<script>
    const trace = {{
        x: {json.dumps(x_pts)},
        y: {json.dumps(y_pts)},
        mode: 'markers',
        type: 'scatter',
        marker: {{ size: 5, color: '#11998e', opacity: 0.6 }}
    }};
    const layout = {{
        xaxis: {{ title: 'Phylogenetic Distance between ASV Pairs' }},
        yaxis: {{ title: 'Niche Distance (' + {json.dumps(column)} + ' Optimum)' }},
        margin: {{ t: 20, b: 50, l: 60, r: 20 }}
    }};
    Plotly.newPlot('mantelScatter', [trace], layout);
</script>
</body>
</html>
"""

    with open(os.path.join(output_dir, "index.html"), "w") as f:
        f.write(index_html)
