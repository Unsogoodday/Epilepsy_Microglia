"""Discover Kumar's sample triplets from GEO SOFT and download MTX/TSV only."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import json
from pathlib import Path
import sys
from urllib.parse import urlparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from epilepsy_microglia.download import fetch, write_json

ACCESSION = 'GSE201048'
SOFT = ACCESSION + '_family.soft.gz'
SOFT_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE201nnn/GSE201048/soft/' + SOFT
ALLOWED = ('.mtx', '.mtx.gz', '.csv', '.csv.gz', '.tsv', '.tsv.gz')


def discover(path):
    with gzip.open(path, 'rt') as f:
        text = f.read()
    if '^SERIES = ' + ACCESSION not in text.splitlines():
        raise ValueError('Unexpected GEO accession')
    series_ids = [l.split(' = ', 1)[1] for l in text.splitlines() if l.startswith('!Series_sample_id = ')]
    sources, samples, filenames = [], set(), set()
    for block in text.split('^SAMPLE = ')[1:]:
        gsm = block.splitlines()[0]
        if gsm in samples:
            raise ValueError('Duplicate GEO sample')
        samples.add(gsm)
        sample_files = []
        for line in block.splitlines():
            if not line.startswith('!Sample_supplementary_file_'):
                continue
            url = line.split(' = ', 1)[1].replace('ftp://', 'https://', 1)
            parsed = urlparse(url)
            if not parsed.path.endswith(ALLOWED):
                continue
            if parsed.scheme != 'https' or parsed.hostname != 'ftp.ncbi.nlm.nih.gov' or f'/{gsm}/suppl/' not in parsed.path:
                raise ValueError('Supplementary URL does not belong to GEO sample')
            name = Path(parsed.path).name
            if name in filenames or not name.startswith(gsm + '_'):
                raise ValueError('Duplicate or conflicting supplementary sample key')
            filenames.add(name)
            sample_files.append(name)
            sources.append(dict(filename=name, url=url, source_sample_id=gsm))
        for suffix in ('_matrix.mtx.gz', '_features.tsv.gz', '_barcodes.tsv.gz'):
            if sum(n.endswith(suffix) for n in sample_files) != 1:
                raise ValueError(f'{gsm}: expected one {suffix}')
        if len(sample_files) != 3:
            raise ValueError(f'{gsm}: unexpected expression file layout')
    if len(samples) != 11 or len(series_ids) != 11 or samples != set(series_ids):
        raise ValueError('Expected the 11 samples listed by GSE201048')
    return sources


def run(root):
    raw, meta = root / 'data/raw/kumar', root / 'data/metadata/kumar'
    staging = meta / 'download_staging'
    raw.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)
    inventory_path = meta / 'download_inventory.json'
    old = {r['filename']: r for r in json.loads(inventory_path.read_text())} if inventory_path.exists() else {}
    item = fetch(dict(filename=SOFT, url=SOFT_URL), raw, staging, old.get(SOFT))
    sources = discover(raw / SOFT)
    allowed = {SOFT, *(s['filename'] for s in sources)}
    completed = {k: v for k, v in old.items() if k in allowed}
    completed[SOFT] = item
    write_json(inventory_path, list(completed.values()))
    with ThreadPoolExecutor(max_workers=3) as pool:
        pending = [pool.submit(fetch, s, raw, staging, old.get(s['filename'])) for s in sources]
        for future in as_completed(pending):
            item = future.result()
            completed[item['filename']] = item
            write_json(inventory_path, sorted(completed.values(), key=lambda r: r['filename']))
    print(f'Kumar: {len(sources)} MTX/TSV files and GEO SOFT verified', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root', type=Path, default=ROOT)
    run(p.parse_args().project_root.resolve())
