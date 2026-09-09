"""Download Ayhan2021/GSE160189 originals; resume, verify, and record provenance.

Python standard library plus curl. No absolute environment-specific paths.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import zipfile
import time
from tempfile import NamedTemporaryFile
from urllib.parse import urljoin

ROOT = Path(__file__).resolve().parents[2]
CLINICAL = '1-s2.0-S0896627321003299-mmc2.xlsx'
# Table S1, linked by the paper associated with GEO PMID 34051145.
# The manifest references sample metadata but does not enumerate this attachment.
CLINICAL_URL = 'https://ars.els-cdn.com/content/image/' + CLINICAL


def write_json(path, value):
    with NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        temporary = Path(f.name)
        json.dump(value, f, indent=2)
        f.write('\n')
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def verify(path):
    if path.name.endswith('.xlsx'):
        with zipfile.ZipFile(path) as z:
            if z.testzip() is not None or 'xl/workbook.xml' not in z.namelist():
                raise ValueError('Invalid Excel workbook')
    if path.name.endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            while f.read(8 * 1024 * 1024):
                pass
    return sha256(path)


def fetch(url, final, staging, expected=None):
    def checked(path):
        # Staging filenames retain the source suffix for format validation.
        digest = verify(path)
        if expected and (digest != expected['sha256'] or path.stat().st_size != expected['size_bytes']):
            raise ValueError(f'Inventory mismatch: {path.name}')
        if final.name.endswith('.html') and 'Index of /geo/series/' not in path.read_text():
            raise ValueError('Not a GEO directory listing')
        return digest
    if final.exists():
        try:
            checked(final)
        except (ValueError, OSError, EOFError, zipfile.BadZipFile):
            final.rename(staging / (final.name + f'.invalid-{time.time_ns()}'))
    if not final.exists():
        partial = staging / final.name
        command = ['curl', '--fail', '--location', '--retry', '5',
                        '--connect-timeout', '30', '--max-time', '1800',
                        '--continue-at', '-', '--silent', '--show-error',
                        '--output', str(partial), url]
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as error:
            if error.returncode not in (33, 36) or not partial.exists():
                raise
            # Server cannot resume: preserve interrupted bytes and restart safely.
            partial.rename(staging / (final.name + f'.unresumable-{time.time_ns()}'))
            subprocess.run(command, check=True)
        try:
            checked(partial)
        except (ValueError, OSError, EOFError, zipfile.BadZipFile):
            partial.rename(staging / (final.name + f'.invalid-{time.time_ns()}'))
            raise ValueError(f'Invalid download preserved in staging: {final.name}; rerun to download afresh')
        os.link(partial, final)
        partial.unlink()
    digest = verify(final)
    print(f'{final.name}: {final.stat().st_size:,} bytes verified', flush=True)
    return dict(filename=final.name, url=url, size_bytes=final.stat().st_size, sha256=digest)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=ROOT)
    p.add_argument('--manifest', type=Path, default=ROOT/'config/dataset_manifest.csv', help='Curated dataset manifest (required; bundled default)')
    a = p.parse_args()
    raw = a.project_root / 'data/raw/ayhan'
    meta = a.project_root / 'data/metadata/ayhan'
    staging = meta / 'download_staging'
    raw.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)
    if a.manifest:
        with a.manifest.open(encoding='utf-8', errors='replace', newline='') as f:
            rows = [r for r in csv.DictReader(f) if r['study_id'] == 'Ayhan2021']
        if len(rows) != 1 or rows[0]['GEO accessions'] != 'GSE160189':
            raise ValueError('Manifest must identify Ayhan2021 as GSE160189')
        write_json(meta / 'manifest_entry.json', rows[0])
    entry = rows[0]
    accession = entry['GEO accessions']
    suppl = entry['GEO supplementary directory']
    base = urljoin(suppl, '../')
    previous = meta/'download_inventory.json'
    old = {r['filename']: r for r in json.loads(previous.read_text())} if previous.exists() else {}
    inventory = [fetch(suppl, raw/'geo_supplementary_index.html', staging, old.get('geo_supplementary_index.html'))]
    names = re.findall(r'href="([^"]+)"', (raw/'geo_supplementary_index.html').read_text())
    names = sorted({n for n in names if n.startswith('GSE160189_') and '/' not in n})
    if 'GSE160189_Hippo_Counts.csv.gz' not in names:
        raise ValueError('Expected matrix missing from GEO listing')
    sources = [(urljoin(suppl, n), n) for n in names]
    sources += [(base+f'soft/{accession}_family.soft.gz', f'{accession}_family.soft.gz'),
                (base+f'matrix/{accession}_series_matrix.txt.gz', f'{accession}_series_matrix.txt.gz'),
                (CLINICAL_URL, CLINICAL)]
    for url, name in sources:
        item = fetch(url, raw/name, staging, old.get(name))
        if name in old and item['sha256'] != old[name]['sha256']:
            raise ValueError(f'Original changed: {name}')
        inventory.append(item)
    write_json(meta/'download_sources.json', inventory)
    write_json(previous, inventory)


if __name__ == '__main__':
    main()
