from __future__ import annotations

import hashlib
import json
import os
import re
import signal
import sqlite3
import subprocess
import sys
import threading
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

def now() -> str:
    return datetime.now(timezone.utc).isoformat()

class JobStore:
    """One process serializes metadata writes; compute workers never write the evidence DB."""
    def __init__(self, database: Path, output_root: Path):
        self.database = Path(database).resolve()
        self.output_root = Path(output_root).resolve()
        self.database.parent.mkdir(parents=True, exist_ok=True)
        self.output_root.mkdir(parents=True, exist_ok=True)
        self.lock = threading.RLock()
        with self.connection() as db:
            db.execute('PRAGMA journal_mode=WAL')
            db.execute('''CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY, algorithm TEXT NOT NULL, action TEXT NOT NULL,
                parameters TEXT NOT NULL, state TEXT NOT NULL, created TEXT NOT NULL,
                updated TEXT NOT NULL, output TEXT NOT NULL, error TEXT,
                code_version TEXT NOT NULL, model_version TEXT, source_manifest TEXT,
                parent_id TEXT)''')
            columns={r[1] for r in db.execute('PRAGMA table_info(jobs)')}
            for name in ('protocol_id','input_hash','output_manifest','output_manifest_sha256'):
                if name not in columns:db.execute('ALTER TABLE jobs ADD COLUMN '+name+' TEXT')

    @contextmanager
    def connection(self):
        # Serialize reads as well as writes in this single metadata process.
        # A read can otherwise race a committing writer and kill the queue thread.
        with self.lock:
            db = sqlite3.connect(self.database, timeout=30)
            db.row_factory = sqlite3.Row
            db.execute('PRAGMA synchronous=NORMAL')
            try:
                with db:
                    yield db
            finally:
                db.close()

    @staticmethod
    def decode(row):
        if row is None:
            raise KeyError('Unknown job')
        result = dict(row)
        result['parameters'] = json.loads(result['parameters'])
        result['source_manifest'] = json.loads(result['source_manifest'] or '[]')
        return result

    def create(self, algorithm, action, parameters, code_version, parent_id=None, protocol_id=None):
        job_id = str(uuid.uuid4())
        output = self.output_root / job_id
        output.mkdir()
        with self.lock, self.connection() as db:
            db.execute('INSERT INTO jobs (id,algorithm,action,parameters,state,created,updated,output,error,code_version,model_version,source_manifest,parent_id,protocol_id) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)', (
                job_id, algorithm, action, json.dumps(parameters), 'queued', now(), now(),
                str(output), None, code_version, None, '[]', parent_id, protocol_id))
        return self.get(job_id)

    def get(self, job_id):
        with self.connection() as db:
            return self.decode(db.execute('SELECT * FROM jobs WHERE id=?', (job_id,)).fetchone())

    def list(self, limit=100):
        with self.connection() as db:
            return [self.decode(r) for r in db.execute('SELECT * FROM jobs ORDER BY created DESC LIMIT ?', (limit,))]

    def pending(self):
        with self.connection() as db:
            return [self.decode(r) for r in db.execute("SELECT * FROM jobs WHERE state='queued' ORDER BY created")]

    def transition(self, job_id, expected, state, error=None):
        with self.lock, self.connection() as db:
            changed = db.execute('UPDATE jobs SET state=?,updated=?,error=? WHERE id=? AND state=?',
                                 (state, now(), error, job_id, expected)).rowcount
        return bool(changed)

    def recover(self):
        # Never silently rerun a partially completed scientific experiment.
        with self.lock, self.connection() as db:
            return db.execute("UPDATE jobs SET state='interrupted',updated=?,error=? WHERE state='running'",
                              (now(), 'Service restarted; explicit retry creates a new run.')).rowcount

    def provenance(self, job_id, model_version, source_manifest, code_version=None, input_hash=None, output_manifest=None, output_manifest_sha256=None):
        with self.lock, self.connection() as db:
            db.execute('UPDATE jobs SET model_version=?,source_manifest=?,code_version=COALESCE(?,code_version),updated=?,input_hash=?,output_manifest=?,output_manifest_sha256=? WHERE id=?',
                       (model_version, json.dumps(source_manifest), code_version, now(), input_hash, output_manifest, output_manifest_sha256, job_id))

    def result(self, job_id):
        job = self.get(job_id)
        if job['state'] != 'succeeded':
            raise ValueError('Result is available only after success')
        files, verification = self.checked_artifacts(job)
        return {'job': job, 'files': files, 'artifact_verification': verification,
                'interpretation': 'Research outputs; predictive scores are not verified treatment effects.'}

    def checked_artifacts(self, job):
        root = Path(job['output']).resolve()
        if not root.is_relative_to(self.output_root):
            raise ValueError('Output escaped configured artifact root')
        files = []
        for path in sorted(root.rglob('*')):
            if not path.is_file() or path.is_symlink():
                continue
            resolved = path.resolve()
            if not resolved.is_relative_to(root):
                continue
            with path.open('rb') as handle:
                digest = hashlib.file_digest(handle, 'sha256').hexdigest()
            files.append({'path': path.relative_to(root).as_posix(), 'bytes': path.stat().st_size, 'sha256': digest})
        manifest_name = job.get('output_manifest')
        if not manifest_name:
            return files, {'frozen_at_success': False, 'status': 'legacy_manifest_unavailable'}
        if not isinstance(manifest_name, str):
            raise ValueError('Invalid output manifest reference')
        manifest_path = (root / manifest_name).resolve()
        if not manifest_path.is_relative_to(root) or not manifest_path.is_file():
            raise ValueError('Output manifest is unavailable inside artifact root')
        manifest_bytes = manifest_path.read_bytes()
        pinned = job.get('output_manifest_sha256')
        if pinned and hashlib.sha256(manifest_bytes).hexdigest() != pinned:
            raise ValueError('Output manifest checksum changed after success')
        manifest = json.loads(manifest_bytes)
        if not isinstance(manifest, dict) or not isinstance(manifest.get('files'), list):
            raise ValueError('Invalid output manifest schema')
        for field in ('protocol_id', 'input_hash', 'model_version'):
            if job.get(field) is not None and manifest.get(field) != job[field]:
                raise ValueError('Output manifest provenance differs: ' + field)
        declared = {}
        for item in manifest['files']:
            if (not isinstance(item, dict) or not isinstance(item.get('path'), str)
                    or not isinstance(item.get('bytes'), int) or isinstance(item['bytes'], bool)
                    or item['bytes'] < 0 or not isinstance(item.get('sha256'), str)):
                raise ValueError('Invalid output manifest artifact record')
            name = item['path']
            if name in declared:
                raise ValueError('Duplicate output manifest artifact')
            declared[name] = {'path': name, 'bytes': item['bytes'], 'sha256': item['sha256']}
        mutable_metadata = {'run.log', manifest_path.relative_to(root).as_posix()}
        actual = {item['path']: item for item in files if item['path'] not in mutable_metadata}
        if actual != declared:
            raise ValueError('Output artifacts differ from frozen manifest')
        return files, {'frozen_at_success': bool(pinned),
                       'status': 'frozen_verified' if pinned else 'legacy_manifest_consistency_only'}

class Runner:
    def __init__(self, store: JobStore, command_builder: Callable):
        self.store = store
        self.command_builder = command_builder
        self.processes: dict[str, subprocess.Popen] = {}
        self.lock = threading.RLock()

    def cancel(self, job_id):
        with self.lock:
            job = self.store.get(job_id)
            if job['state'] in ('succeeded', 'failed', 'cancelled', 'interrupted'):
                return job
            self.store.transition(job_id, job['state'], 'cancelled')
            process = self.processes.get(job_id)
            if process and process.poll() is None:
                if os.name == 'posix':
                    os.killpg(process.pid, signal.SIGTERM)
                else:
                    process.terminate()
        return self.store.get(job_id)

    def run(self, job_id):
        with self.lock:
            if not self.store.transition(job_id, 'queued', 'running'):
                return self.store.get(job_id)
            job = self.store.get(job_id)
            output = Path(job['output'])
            log = None
            try:
                command = self.command_builder(job)
                if not command or not all(isinstance(x, str) for x in command):
                    raise ValueError('Invalid command')
                (output / 'command.json').write_text(json.dumps(command, indent=2), encoding='utf-8')
                log = (output / 'run.log').open('w', encoding='utf-8')
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                    stdin=subprocess.DEVNULL, start_new_session=os.name == 'posix', shell=False)
                self.processes[job_id] = process
            except Exception as exc:
                if log is not None:
                    log.close()
                self.store.transition(job_id, 'running', 'failed', str(exc))
                return self.store.get(job_id)
        code = process.wait()
        log.close()
        with self.lock:
            self.processes.pop(job_id, None)
            if code != 0:
                self.store.transition(job_id, 'running', 'failed', f'Process exited with code {code}')
                return self.store.get(job_id)
            provenance = output / 'provenance.json'
            if job.get('protocol_id') and not provenance.is_file():
                self.store.transition(job_id, 'running', 'failed', 'Frozen protocol output provenance is missing')
                return self.store.get(job_id)
            if provenance.exists():
                try:
                    value = json.loads(provenance.read_text(encoding='utf-8'))
                    if not isinstance(value, dict) or not isinstance(value.get('sources', []), list):
                        raise ValueError('Invalid result provenance schema')
                    for field in ('model_version', 'code_version', 'input_hash', 'protocol_id'):
                        if value.get(field) is not None and not isinstance(value[field], str):
                            raise ValueError('Invalid result provenance field: ' + field)
                    manifest_name = value.get('output_manifest')
                    if job.get('protocol_id') and (not manifest_name or not value.get('input_hash')):
                        raise ValueError('Frozen protocol output manifest and input hash are required')
                    if job.get('protocol_id'):
                        if not re.fullmatch(r'[a-f0-9]{64}', value['input_hash']):
                            raise ValueError('Frozen protocol input hash must be SHA256')
                        if value.get('protocol_id') != job['protocol_id']:
                            raise ValueError('Result provenance protocol differs from queued job')
                    manifest_hash = None
                    if manifest_name is not None:
                        if not isinstance(manifest_name, str):
                            raise ValueError('Invalid output manifest reference')
                        manifest_path = (output / manifest_name).resolve()
                        if not manifest_path.is_relative_to(output.resolve()):
                            raise ValueError('Output manifest escaped artifact root')
                        manifest_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
                    self.store.provenance(job_id, value.get('model_version'), value.get('sources', []),value.get('code_version'),
                                          value.get('input_hash'),manifest_name,manifest_hash)
                    self.store.checked_artifacts(self.store.get(job_id))
                except (OSError, ValueError, TypeError, KeyError) as exc:
                    self.store.transition(job_id, 'running', 'failed', 'Invalid result provenance: ' + str(exc))
                    return self.store.get(job_id)
            self.store.transition(job_id, 'running', 'succeeded')
        return self.store.get(job_id)
