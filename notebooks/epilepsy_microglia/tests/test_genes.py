import tempfile
import unittest
import anndata as ad
import numpy as np
from test_contract import example
from epilepsy_microglia.genes import standardize_var_names
from epilepsy_microglia.io import write_study
from epilepsy_microglia.validation import validate_adata


class GeneTests(unittest.TestCase):
    def test_writer_converts_and_preserves_counts(self):
        a = example()
        a.var_names = ['ENSG000001.2', 'ENSG000002']
        a.var['source_gene_id'] = a.var_names.copy()
        a.var['gene_symbol'] = ['TBCE', 'TBCE']
        self.assertTrue(validate_adata(a, study='kumar'))
        with tempfile.TemporaryDirectory() as root:
            result = ad.read_h5ad(write_study(a, study='kumar', project_root=root))
        self.assertEqual(list(result.var_names), ['TBCE', 'TBCE-1'])
        self.assertEqual(list(result.var.source_gene_id), ['ENSG000001.2', 'ENSG000002'])
        np.testing.assert_array_equal(result.X.toarray(), a.X.toarray())
        np.testing.assert_array_equal(result.layers['counts'].toarray(), a.X.toarray())
        self.assertEqual(validate_adata(result, study='kumar'), [])
        standardize_var_names(result)
        self.assertEqual(list(result.var_names), ['TBCE', 'TBCE-1'])

    def test_missing_mapping_fails_without_mutation(self):
        for symbol in [None, '', ' ', 'ENSG000001']:
            a = example()
            a.var_names = ['ENSG000001', 'existing']
            a.var['gene_symbol'] = [symbol, None]
            with self.assertRaises(ValueError):
                standardize_var_names(a)
            self.assertEqual(list(a.var_names), ['ENSG000001', 'existing'])

    def test_mixed_names_and_source_preservation(self):
        a = example()
        a.var_names = ['ENSG000001', 'existing']
        del a.var['source_gene_id']
        a.var['gene_symbol'] = ['GENE', None]
        standardize_var_names(a)
        self.assertEqual(list(a.var_names), ['GENE', 'existing'])
        self.assertEqual(list(a.var.source_gene_id), ['ENSG000001', 'existing'])
