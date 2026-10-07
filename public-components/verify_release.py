"""Verify release files without network, model execution or private inputs."""
import hashlib,json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'PUBLIC_MANIFEST.json').read_text())
for item in manifest['files']:
    relative=Path(item['path'])
    if relative.is_absolute() or '..' in relative.parts:raise ValueError('Invalid manifest path')
    raw=(ROOT/relative).read_bytes()
    if len(raw)!=item['bytes'] or hashlib.sha256(raw).hexdigest()!=item['sha256']:
        raise ValueError('Distribution file differs: '+item['path'])
source=ROOT/'SOURCE_MANIFEST.json'
provenance=json.loads((ROOT/'public-api/provenance.json').read_text())
if hashlib.sha256(source.read_bytes()).hexdigest()!=provenance['source_manifest_sha256']:
    raise ValueError('API source manifest pin differs')
for item in json.loads(source.read_text())['files']:
    if hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()!=item['distribution_sha256']:
        raise ValueError('Source distribution differs: '+item['path'])
print(json.dumps({'status':'passed','distribution_files_verified':len(manifest['files']),
                  'source_manifest_sha256':provenance['source_manifest_sha256']}))
