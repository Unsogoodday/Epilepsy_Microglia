import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from epilepsy_microglia.metadata import standardize_obs, join_metadata
from epilepsy_microglia.validation import validate_adata
from epilepsy_microglia.io import write_study


def example(sparse_counts=True):
    x = np.array([[0, 2], [3, 0]], dtype=np.int32)
    if sparse_counts:
        x = sparse.csr_matrix(x)
    obs = standardize_obs(pd.DataFrame(index=['cell1', 'cell2']), 'kumar')
    obs['source_file'] = 'kumar/original.mtx.gz'
    a = ad.AnnData(X=x, obs=obs, var=pd.DataFrame({'source_gene_id': ['g1', 'g2']}, index=['g1', 'g2']))
    a.layers['counts'] = x.copy()
    a.uns['ingestion'] = dict(raw_counts_available=True, expression_type='raw_counts', source_files=['kumar/original.mtx.gz'], metadata_sources=[])
    return a


class ContractTests(unittest.TestCase):
    def test_counts_roundtrip_and_no_overwrite(self):
        for use_sparse in [False, True]:
            with self.subTest(sparse=use_sparse), tempfile.TemporaryDirectory() as d:
                a = example(use_sparse)
                self.assertEqual(validate_adata(a, study='kumar'), [])
                path = write_study(a, study='kumar', project_root=d)
                result = ad.read_h5ad(path)
                self.assertEqual(validate_adata(result, study='kumar'), [])
                self.assertTrue(result.obs.fcd_subtype.isna().all())
                self.assertEqual(list(result.obs.source_cell_id), ['cell1', 'cell2'])
                original = path.read_bytes()
                with self.assertRaises(FileExistsError):
                    write_study(a, study='kumar', project_root=d)
                self.assertEqual(path.read_bytes(), original)
                self.assertFalse((Path(d)/'data/raw').exists())

    def test_missing_clinical_and_explicit_join(self):
        obs = pd.DataFrame({'source_sample_id': ['b', 'a', 'z']}, index=['c3', 'c1', 'c2'])
        lookup = pd.DataFrame({'source_sample_id': ['a', 'b'], 'diagnosis': ['FCD', 'TLE']})
        joined = join_metadata(obs, lookup, key='source_sample_id')
        self.assertEqual(list(joined.index), list(obs.index))
        self.assertTrue(pd.isna(joined.loc['c2', 'diagnosis']))
        standardized = standardize_obs(joined, 'kumar')
        self.assertTrue(standardized.fcd_subtype.isna().all())
        self.assertEqual(standardized.loc['c1', 'diagnosis'], 'FCD')
        with self.assertRaises(ValueError):
            join_metadata(obs, pd.concat([lookup, lookup]), key='source_sample_id')
        with self.assertRaises(ValueError):
            join_metadata(joined, lookup, key='source_sample_id')

    def test_bad_counts(self):
        for value in [-1, 0.5, float('nan'), float('inf')]:
            a = example(False)
            a.layers['counts'] = np.array([[value, 2], [3, 0]])
            self.assertTrue(validate_adata(a, study='kumar'))
        a = example()
        a.X = a.X * 2
        self.assertTrue(validate_adata(a, study='kumar'))

    def test_missing_counts_and_source_ids(self):
        a = example()
        del a.layers['counts']
        self.assertTrue(validate_adata(a, study='kumar'))
        a.uns['ingestion'].update(raw_counts_available=False, expression_type='source_normalized')
        self.assertEqual(validate_adata(a, study='kumar'), [])
        del a.var['source_gene_id']
        self.assertTrue(validate_adata(a, study='kumar'))

    def test_duplicate_ids_and_study(self):
        a = example()
        a.obs_names = ['same', 'same']
        self.assertTrue(validate_adata(a, study='kumar'))
        self.assertTrue(validate_adata(example(), study='liu'))
        with self.assertRaises(ValueError):
            standardize_obs(example().obs, 'liu')

    def test_symlink_refused(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'data').mkdir()
            (root/'data/raw').mkdir()
            (root/'data/processed').symlink_to(root/'data/raw', target_is_directory=True)
            with self.assertRaises(ValueError):
                write_study(example(), study='kumar', project_root=root)
            self.assertEqual(list((root/'data/raw').iterdir()), [])


if __name__ == '__main__':
    unittest.main()
