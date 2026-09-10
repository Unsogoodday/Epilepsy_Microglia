"""Download Kumar MTX/TSV triplets only; run offline ingestion separately.

The user's MTX/tabular-only scope supersedes the manifest's all-supplements rule.
Completed raw files are verified in place and never overwritten.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
from epilepsy_microglia.download import fetch, write_json


def source_list(root):
    with (root / 'config/dataset_manifest.csv').open(errors='replace') as f:
        entries = [r for r in csv.DictReader(f) if r['study_id'] == 'Kumar2022']
    if len(entries) != 1 or entries[0]['GEO accessions'] != 'GSE201048':
        raise ValueError('Expected one Kumar2022/GSE201048 manifest entry')
    with (root / 'data/metadata/kumar/source_triplets.csv').open() as f:
        sources = list(csv.DictReader(f))
    seen = set()
    for source in sources:
        name, gsm = source['filename'], source['source_sample_id']
        if (not re.fullmatch(r'GSM60496(?:3[2-9]|4[0-2])_Sample\d+_(?:matrix\.mtx|features\.tsv|barcodes\.tsv)\.gz', name)
                or not name.startswith(gsm + '_') or name in seen):
            raise ValueError(f'Unexpected or duplicate triplet source: {name}')
        seen.add(name)
        expected_url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6049nnn/{gsm}/suppl/{name}'
        if source['url'] != expected_url:
            raise ValueError(f'Unexpected source URL: {name}')
        source['expected_size_bytes'] = int(source['expected_size_bytes'])
    if len(sources) != 33 or len({s['source_sample_id'] for s in sources}) != 11:
        raise ValueError('Expected exactly 11 MTX triplets')
    return sources


def run(root):
    sources = source_list(root)
    raw, meta = root / 'data/raw/kumar', root / 'data/metadata/kumar'
    staging = meta / 'download_staging'
    raw.mkdir(parents=True, exist_ok=True)
    staging.mkdir(parents=True, exist_ok=True)
    previous = meta / 'download_inventory.json'
    allowed = {s['filename'] for s in sources}
    old = {r['filename']: r for r in json.loads(previous.read_text()) if r['filename'] in allowed} if previous.exists() else {}
    with ThreadPoolExecutor(max_workers=3) as pool:
        pending = [pool.submit(fetch, source, raw, staging, old.get(source['filename'])) for source in sources]
        for future in as_completed(pending):
            item = future.result()
            old[item['filename']] = item
            write_json(previous, sorted(old.values(), key=lambda r: r['filename']))
    print(f'Complete: {len(sources)} MTX/TSV files verified', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=ROOT)
    run(parser.parse_args().project_root.resolve())
