"""Immutable file identities and source registry hooks; no implicit downloads."""
from __future__ import annotations
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(2 ** 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def registry_hook(path: Path | None, record: dict) -> None:
    """Append/update by source ID; unchanged inputs cannot silently change identity."""
    if path is None:
        return
    registry = json.loads(path.read_text(encoding='utf-8')) if path.exists() else {'schema_version': 1, 'sources': []}
    if not isinstance(registry.get('sources'), list):
        raise ValueError('registry must contain a sources list')
    existing = next((r for r in registry['sources'] if r['id'] == record['id']), None)
    if existing is not None and existing.get('sha256') != record.get('sha256'):
        raise ValueError('conflicting registry source identity')
    if existing is None:
        registry['sources'].append(record)
    else:
        existing.update(record)
    atomic_json(path, registry)


def file_record(path: Path, source: str, url: str, version: str | None, rights: str) -> dict:
    digest = sha256(path)
    return {'id': f'{source}:{digest}', 'source': source, 'url': url,
            'path': str(path.resolve()), 'sha256': digest, 'bytes': path.stat().st_size,
            'version': version, 'recorded_at': now(), 'rights': rights,
            'access': 'local_snapshot', 'redistribution_verified': False}
