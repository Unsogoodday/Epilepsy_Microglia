"""Download only the inspected Galvao source list; resume in staging, verify gzip.

Run separately from ingestion. Existing raw files are never changed.
"""
from pathlib import Path
import concurrent.futures
import gzip
import hashlib
import json
import os
import subprocess

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / 'data/raw/galvao'
META = ROOT / 'data/metadata/galvao'
STAGING = META / 'download_staging'


def verify(path):
    if path.suffix == '.gz':
        with gzip.open(path, 'rb') as f:
            while f.read(8 * 1024 * 1024):
                pass  # Read full stream: gzip CRC and length validation.
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def download(source):
    name = source['filename']
    if Path(name).name != name:
        raise ValueError(name)
    final = RAW / name
    if not final.exists():
        partial = STAGING / name
        subprocess.run(['curl', '--fail', '--location', '--retry', '5',
                        '--retry-delay', '3', '--connect-timeout', '30',
                        '--continue-at', '-', '--silent', '--show-error',
                        '--output', str(partial), source['url']], check=True)
        if "expected_size_bytes" in source and partial.stat().st_size != source["expected_size_bytes"]:
            raise ValueError(f"Unexpected size for {name}")
        digest = verify(partial)
        os.link(partial, final)  # Fails instead of overwriting an original.
        partial.unlink()
    else:
        digest = verify(final)
    result = dict(source, size_bytes=final.stat().st_size, sha256=digest,
                  integrity='full gzip stream CRC/length passed' if name.endswith('.gz') else 'SHA256 recorded')
    print(f'{name}: {result["size_bytes"]:,} bytes, verified', flush=True)
    return result


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    STAGING.mkdir(parents=True, exist_ok=True)
    sources = json.loads((META / 'download_sources.json').read_text())
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(download, sources))
    (META / 'download_inventory.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
