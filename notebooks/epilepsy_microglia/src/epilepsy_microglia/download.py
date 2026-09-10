"""Immutable raw downloads with resumable staging and recorded integrity."""
import gzip
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def verify(path):
    if path.name.endswith('.gz'):
        with gzip.open(path, 'rb') as f:
            while f.read(8 * 1024 * 1024):
                pass
    return sha256(path)


def fetch(source, raw, staging, expected=None):
    name = source['filename']
    if Path(name).name != name:
        raise ValueError(name)
    final, partial = raw / name, staging / name

    def checked(path):
        digest = verify(path)
        size = source.get('expected_size_bytes')
        if size is not None and path.stat().st_size != size:
            raise ValueError(f'Size mismatch: {name}')
        if expected and (digest != expected['sha256'] or path.stat().st_size != expected['size_bytes']):
            raise ValueError(f'Inventory mismatch: {name}')
        return digest

    if not final.exists():
        command = ['curl', '--fail', '--location', '--retry', '5', '--connect-timeout', '30',
                   '--continue-at', '-', '--silent', '--show-error', '--output', str(partial), source['url']]
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            if e.returncode not in (33, 36) or not partial.exists():
                raise
            partial.rename(staging / (name + f'.unresumable-{time.time_ns()}'))
            subprocess.run(command, check=True)
        try:
            digest = checked(partial)
        except (ValueError, OSError, EOFError):
            partial.rename(staging / (name + f'.invalid-{time.time_ns()}'))
            raise
        os.link(partial, final)
        partial.unlink()
    else:
        digest = checked(final)  # Never move, replace, or repair completed raw files.
    print(f'{name}: {final.stat().st_size:,} bytes verified', flush=True)
    return dict(source, size_bytes=final.stat().st_size, sha256=digest)
