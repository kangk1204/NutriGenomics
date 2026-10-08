"""Verify release files without network, model execution or private inputs."""
import argparse,hashlib,json,subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parent
manifest=json.loads((ROOT/'PUBLIC_MANIFEST.json').read_text())
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--tracked',action='store_true',help='Also require exact coverage of Git-tracked public-components files.')
args=parser.parse_args()
names=[item['path'] for item in manifest['files']]
if len(names)!=len(set(names)):raise ValueError('Duplicate distribution manifest path')
if args.tracked:
    git_root=Path(subprocess.check_output(['git','rev-parse','--show-toplevel'],cwd=ROOT,text=True).strip()).resolve()
    prefix=ROOT.relative_to(git_root).as_posix()+'/'
    tracked=subprocess.check_output(['git','ls-files','-z','--',prefix],cwd=git_root).decode('utf-8').split('\0')
    expected={name.removeprefix(prefix) for name in tracked if name}
    if set(names)|{'PUBLIC_MANIFEST.json'}!=expected:
        raise ValueError('Distribution manifest does not cover exactly the Git-tracked release files')
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
