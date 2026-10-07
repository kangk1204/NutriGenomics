import hashlib
import json
from pathlib import Path
import pytest
from fastapi.testclient import TestClient
from nutriomics_final.api import create_app
from nutriomics_final.research_registry import ResearchRegistry


def bundle(tmp_path):
    folder = tmp_path / 'snapshot'
    folder.mkdir()
    payload = b'{"official_target":null,"observed_auroc":0.925,"official_achievement_percent":null}'
    (folder / 'goals.json').write_bytes(payload)
    manifest = {"schema_version": 1, "snapshot_version": "test", "generated_utc": "2026-10-03T00:00:00Z", "artifacts": {
        "goals": {"path": "goals.json", "sha256": hashlib.sha256(payload).hexdigest()}}}
    path = folder / 'manifest.json'
    path.write_text(json.dumps(manifest))
    return path


def test_status_preserves_missing_official_target(tmp_path):
    path = bundle(tmp_path)
    with TestClient(create_app(tmp_path, research_manifest=path)) as client:
        result = client.get('/research/status')
        assert result.status_code == 200
        assert result.json()['goals']['official_achievement_percent'] is None
        assert client.get('/research/evidence/unknown').status_code == 404


def test_changed_evidence_fails_closed(tmp_path):
    path = bundle(tmp_path)
    (path.parent / 'goals.json').write_text('{}')
    with TestClient(create_app(tmp_path, research_manifest=path)) as client:
        assert client.get('/research/evidence/goals').status_code == 409


def test_manifest_and_evidence_path_escape_rejected(tmp_path):
    with pytest.raises(ValueError, match='manifest'):
        ResearchRegistry(tmp_path, tmp_path.parent / 'outside.json')
    path = bundle(tmp_path)
    manifest = json.loads(path.read_text())
    manifest['artifacts']['goals']['path'] = '../outside.json'
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='escaped'):
        ResearchRegistry(tmp_path, path).evidence('goals')


def test_snapshot_missing_has_explicit_service_status(tmp_path):
    with TestClient(create_app(tmp_path)) as client:
        assert client.get('/research/status').status_code == 503
