"""Download allowed Ayhan expression files discovered in original GEO SOFT.

Only MTX/CSV/TSV expression and GEO SOFT metadata; no clinical spreadsheets.
"""
import argparse
import gzip
import json
from pathlib import Path
import sys
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from epilepsy_microglia.download import fetch, write_json

ACCESSION = 'GSE160189'
SOFT = ACCESSION + '_family.soft.gz'
SOFT_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160189/soft/' + SOFT
ALLOWED = ('.mtx', '.mtx.gz', '.csv', '.csv.gz', '.tsv', '.tsv.gz')


def discover(path):
    with gzip.open(path, 'rt') as f:
        lines = f.read().splitlines()
    if '^SERIES = ' + ACCESSION not in lines:
        raise ValueError('Unexpected GEO accession')
    sources = {}
    for line in lines:
        if line.startswith(('!Series_supplementary_file = ', '!Sample_supplementary_file_')):
            url = line.split(' = ', 1)[1].replace('ftp://', 'https://', 1)
            parsed = urlparse(url)
            if not parsed.path.endswith(ALLOWED):
                continue
            if parsed.scheme != 'https' or parsed.hostname != 'ftp.ncbi.nlm.nih.gov':
                raise ValueError('Unexpected GEO supplementary host')
            name = Path(parsed.path).name
            source = dict(filename=name, url=url)
            if name in sources and sources[name] != source:
                raise ValueError('Duplicate source filename')
            sources[name] = source
    if not sources:
        raise ValueError('No allowed expression sources in GEO')
    return list(sources.values())


def run(root):
    raw, meta = root / 'data/raw/ayhan', root / 'data/metadata/ayhan'
    staging = meta / 'download_staging'
    raw.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)
    inventory_path = meta / 'download_inventory.json'
    old = {r['filename']: r for r in json.loads(inventory_path.read_text())} if inventory_path.exists() else {}
    soft_source = dict(filename=SOFT, url=SOFT_URL)
    inventory = [fetch(soft_source, raw, staging, old.get(SOFT))]
    sources = discover(raw / SOFT)
    # The SOFT discovery document is the source authority, not this audit output.
    allowed = {SOFT, *(s['filename'] for s in sources)}
    completed = {k: v for k, v in old.items() if k in allowed}
    completed[SOFT] = inventory[0]
    write_json(inventory_path, list(completed.values()))
    for source in sources:
        item = fetch(source, raw, staging, old.get(source['filename']))
        inventory.append(item)
        completed[item['filename']] = item
        write_json(inventory_path, list(completed.values()))
    write_json(inventory_path, inventory)
    print(f'Ayhan: {len(sources)} expression file(s) and GEO SOFT verified', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=ROOT)
    run(p.parse_args().project_root.resolve())
