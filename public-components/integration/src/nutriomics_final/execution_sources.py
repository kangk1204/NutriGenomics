"""Bind first-party subprocess imports to one auditable source layout."""
import hashlib
import json
from pathlib import Path
import sys

PACKAGES = {'evidence':'nutriomics-evidence-engine', 'methylation':'nutriomics-htn-methylation',
            'intervention':'nutriomics-intervention-atlas', 'compound-target':'nutriomics-compound-target'}
MODULES = {'evidence':'nutriomics_evidence', 'methylation':'nutriomics_methylation',
           'intervention':'nutriomics_atlas', 'compound-target':'nutriomics_dti'}
PUBLIC = {'evidence':'evidence', 'methylation':'methylation',
          'intervention':'intervention-atlas', 'compound-target':'compound-target'}


def source_layout(root, algorithm):
    root=Path(root).resolve()
    definitions=(('nutriomics-final-report',PACKAGES[algorithm]),('integration',PUBLIC[algorithm]))
    modules=('nutriomics_final',MODULES[algorithm])
    matches=[]
    for names in definitions:
        entries=[]
        for name,module in zip(names,modules,strict=True):
            directory=(root/name).resolve()
            src=(directory/'src').resolve()
            package=(src/module).resolve()
            if not directory.is_relative_to(root) or not src.is_relative_to(root) or not package.is_relative_to(root):
                raise ValueError('Auditable source escaped research root: '+name)
            if not (package/'__init__.py').is_file():break
            entries.append({'package':'nutriomics-final-report' if module=='nutriomics_final' else PACKAGES[algorithm],
                            'directory':directory,'src':src,'module':module})
        if len(entries)==2:matches.append(entries)
    if not matches:raise ValueError('Auditable source checkout missing for '+algorithm)
    if len(matches)!=1:raise ValueError('Ambiguous auditable source layout for '+algorithm)
    return matches[0]


_BOOTSTRAP = r"""import hashlib,json,runpy,sys
from pathlib import Path
configuration=json.loads(sys.argv[1]);module=sys.argv[2];arguments=sys.argv[3:]
roots={name:Path(path).resolve() for name,path in configuration['roots'].items()}
sys.path[:0]=[str(path.parent) for path in roots.values()]
sys.argv=[module,*arguments]
namespace=runpy.run_module(module,run_name='__main__',alter_sys=True)
origins={}
for name,value in list(sys.modules.items()):
    prefix=name.split('.')[0]
    if prefix not in roots:continue
    origin=getattr(value,'__file__',None)
    if not origin:raise ValueError('Unverifiable first-party module: '+name)
    path=Path(origin).resolve()
    if not path.is_relative_to(roots[prefix]) or path.suffix!='.py':raise ValueError('First-party module escaped source snapshot: '+name)
    origins[name]={'module':name,'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
path=Path(namespace['__file__']).resolve();prefix=module.split('.')[0]
if not path.is_relative_to(roots[prefix]) or path.suffix!='.py':raise ValueError('Entrypoint escaped source snapshot')
origins[module]={'module':module,'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
if configuration['receipt']:
    Path(configuration['receipt']).write_text(json.dumps({'module':module,'origins':[origins[name] for name in sorted(origins)]},indent=2),encoding='utf-8')
"""


def bound_command(root, algorithm, command, receipt=None):
    if len(command)<3 or command[1]!='-m':raise ValueError('Expected an allowlisted Python module command')
    roots={entry['module']:str(entry['src']/entry['module']) for entry in source_layout(root,algorithm)}
    if command[2].split('.')[0] not in roots:raise ValueError('Unexpected first-party entrypoint')
    configuration={'roots':roots,'receipt':str(receipt) if receipt is not None else None}
    # Isolated mode ignores inherited PYTHONPATH, the current directory and user site.
    return [command[0],'-I','-c',_BOOTSTRAP,json.dumps(configuration),*command[2:]]


def current_module_origins(root,algorithm):
    roots={entry['module']:(entry['src']/entry['module']).resolve() for entry in source_layout(root,algorithm)}
    values={}
    for name,module in list(sys.modules.items()):
        prefix=name.split('.')[0]
        if prefix not in roots:continue
        origin=getattr(module,'__file__',None)
        if not origin:raise ValueError('Unverifiable first-party module: '+name)
        path=Path(origin).resolve()
        if not path.is_relative_to(roots[prefix]) or path.suffix!='.py':raise ValueError('First-party module escaped source snapshot: '+name)
        values[name]={'module':name,'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    return [values[name] for name in sorted(values)]


def verify_module_origins(snapshots,origins):
    declared={str(Path(item['path']).resolve()):item['sha256'] for snapshot in snapshots for item in snapshot['files']}
    if not origins:raise ValueError('Executed first-party module receipt is empty')
    for item in origins:
        if declared.get(str(Path(item['path']).resolve()))!=item['sha256']:
            raise ValueError('Executed module differs from audited source snapshot: '+item['module'])
