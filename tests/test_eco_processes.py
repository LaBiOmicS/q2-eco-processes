import os
import tempfile
import unittest
import numpy as np
import pandas as pd
import biom
import skbio
import qiime2
from io import StringIO
from q2_eco_processes._processes import calculate_processes, calculate_matrices
from q2_eco_processes._visualizer import summarize_processes
from q2_eco_processes._signal import check_phylogenetic_signal

class TestEcoProcessesComprehensive(unittest.TestCase):

    def setUp(self):
        # 10 ASVs x 6 samples across 2 groups (Control vs Treated)
        np.random.seed(42)
        matrix = np.random.randint(0, 500, size=(10, 6))
        matrix[0, :2] = 0

        self.sample_ids = [f'Sample{i+1}' for i in range(6)]
        self.asv_ids = [f'ASV{i+1}' for i in range(10)]

        self.table = biom.Table(
            matrix,
            observation_ids=self.asv_ids,
            sample_ids=self.sample_ids
        )

        tree_nwk = "(((ASV1:0.1,ASV2:0.2)n1:0.1,(ASV3:0.15,ASV4:0.25)n2:0.2)n3:0.1,((ASV5:0.1,ASV6:0.2)n4:0.1,(ASV7:0.3,(ASV8:0.1,(ASV9:0.05,ASV10:0.05)n5:0.1)n6:0.1)n7:0.2)n8:0.3);"
        self.tree = skbio.TreeNode.read(StringIO(tree_nwk))

        meta_df = pd.DataFrame({
            'Group': ['Control', 'Control', 'Control', 'Treated', 'Treated', 'Treated'],
            'pH': [6.5, 6.7, 6.8, 8.1, 8.3, 8.5]
        }, index=self.sample_ids)
        meta_df.index.name = 'sample-id'
        self.metadata = qiime2.Metadata(meta_df)

    def test_missing_tree_raises(self):
        with self.assertRaises(ValueError):
            calculate_processes(self.table, tree=None)

    def test_empty_table_raises(self):
        empty_table = biom.Table(np.zeros((10, 6)), observation_ids=self.asv_ids, sample_ids=self.sample_ids)
        with self.assertRaises(ValueError):
            calculate_processes(empty_table, tree=self.tree)

    def test_mismatched_asvs_raises(self):
        tree_nwk = "(ASV99:0.1,ASV100:0.2);"
        tree = skbio.TreeNode.read(StringIO(tree_nwk))
        with self.assertRaises(ValueError):
            calculate_processes(self.table, tree=tree)

    def test_sample_level_calculation(self):
        proc_table = calculate_processes(
            table=self.table,
            tree=self.tree,
            permutations=49,
            min_frequency=10,
            n_jobs=1
        )
        self.assertIsInstance(proc_table, biom.Table)
        self.assertEqual(len(proc_table.ids(axis='sample')), 6)

    def test_calculate_matrices(self):
        beta_nti_dm, rc_bray_dm = calculate_matrices(
            table=self.table,
            tree=self.tree,
            permutations=19,
            n_jobs=1
        )
        self.assertIsInstance(beta_nti_dm, skbio.DistanceMatrix)
        self.assertIsInstance(rc_bray_dm, skbio.DistanceMatrix)
        self.assertEqual(beta_nti_dm.shape, (6, 6))
        self.assertEqual(rc_bray_dm.shape, (6, 6))

    def test_check_phylogenetic_signal(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            check_phylogenetic_signal(
                output_dir=tmpdir,
                table=self.table,
                tree=self.tree,
                metadata=self.metadata,
                column='pH',
                permutations=49
            )
            index_path = os.path.join(tmpdir, "index.html")
            self.assertTrue(os.path.exists(index_path))
            with open(index_path, 'r') as f:
                content = f.read()
                self.assertIn("Phylogenetic Signal Analysis", content)
                self.assertIn("Mantel Statistic", content)

    def test_visualizer(self):
        proc_table = calculate_processes(
            table=self.table,
            tree=self.tree,
            permutations=19,
            n_jobs=1
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            summarize_processes(tmpdir, proc_table)
            self.assertTrue(os.path.exists(os.path.join(tmpdir, "index.html")))

if __name__ == '__main__':
    unittest.main()
