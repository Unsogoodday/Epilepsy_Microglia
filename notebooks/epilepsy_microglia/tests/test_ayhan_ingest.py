"""Small adversarial checks for the study-specific count parser and metadata join."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import hashlib

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('ayhan_ingest', ROOT/'scripts/ingest/ayhan.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class AyhanTests(unittest.TestCase):
    def test_download_reuses_verified_file_and_repairs_corruption(self):
        import ayhan as download
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            staging = root/'staging'
            staging.mkdir()
            final = root/'counts.csv.gz'
            payload = gzip.compress(b'gene,A56_A\nG,1\n')
            final.write_bytes(payload)
            expected = dict(size_bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
            with patch.object(download.subprocess, 'run') as call:
                download.fetch('https://example.test/counts.csv.gz', final, staging, expected)
                call.assert_not_called()
            final.write_bytes(b'<html>download failed</html>')
            def complete(command, check):
                Path(command[command.index('--output')+1]).write_bytes(payload)
            with patch.object(download.subprocess, 'run', side_effect=complete):
                download.fetch('https://example.test/counts.csv.gz', final, staging, expected)
            self.assertEqual(final.read_bytes(), payload)
            self.assertEqual(len(list(staging.glob('*.invalid-*'))), 1)

    def test_download_rejects_truncated_gzip_before_publication(self):
        import ayhan as download
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            staging = root/'staging'
            staging.mkdir()
            def truncated(command, check):
                Path(command[command.index('--output')+1]).write_bytes(gzip.compress(b'abc')[:-5])
            with patch.object(download.subprocess, 'run', side_effect=truncated), self.assertRaises(ValueError):
                download.fetch('https://example.test/counts.csv.gz', root/'counts.csv.gz', staging)
            self.assertFalse((root/'counts.csv.gz').exists())
            self.assertEqual(len(list(staging.glob('*.invalid-*'))), 1)

    def parse(self, text):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'counts.csv.gz'
            with gzip.open(p, 'wt') as f:
                f.write(text)
            return m.read_counts(p)

    def test_orientation_zero_genes_and_cell_order(self):
        x, cells, genes = self.parse('gene,P57_B,A56_A\nG1,0,3\nG2,2,0\nG3,0,0\n')
        np.testing.assert_array_equal(x.toarray(), [[0,2,0],[3,0,0]])
        self.assertEqual(list(cells), ['P57_B','A56_A'])
        self.assertEqual(genes, ['G1','G2','G3'])

    def test_invalid_counts_rejected(self):
        for values in ['1,-1','0,1.5','0,NaN','0,','1,2,3','1,2147483648']:
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.parse('gene,A56_A,P57_B\nG1,'+values+'\n')

    def test_duplicate_axes_rejected(self):
        for text in ['gene,A56_A,A56_A\nG,1,2\n', 'gene,A56_A\nG,1\nG,2\n']:
            with self.assertRaises(ValueError):
                self.parse(text)

    def test_geo_and_clinical_join(self):
        raw = ROOT/'data/raw/ayhan'
        if not (raw/m.CLINICAL).exists():
            self.skipTest('Run downloader to enable real-source metadata checks')
        samples = m.read_clinical(raw/m.CLINICAL, m.read_geo(raw/'GSE160189_family.soft.gz'))
        self.assertEqual(len(samples), 10)
        self.assertEqual(samples.donor_id.nunique(), 5)
        self.assertEqual(samples.loc['A56','donor_id'], 'Donor1')
        self.assertEqual(samples.loc['P57','seizure_frequency_per_month'], 3)
        self.assertEqual(set(samples.index[samples.hemisphere_report_conflict]), {'A76','P76'})
        self.assertTrue(samples.medications.notna().all())


if __name__ == '__main__':
    unittest.main()
