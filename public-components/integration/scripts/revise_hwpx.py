"""Modify exact paragraph text in a copied HWPX; preserve other package members."""
import argparse
from collections import Counter
import io
import json
import os
import tempfile
from pathlib import Path
import zipfile
from lxml import etree as ET

HP='http://www.hancom.co.kr/hwpml/2011/paragraph'

def own_text_nodes(paragraph):
    return paragraph.findall(f'./{{{HP}}}run/{{{HP}}}t')

def own_text(paragraph):
    return ''.join(''.join(t.itertext()) for t in own_text_nodes(paragraph)).strip()

def paragraphs(path):
    rows=[]
    with zipfile.ZipFile(path) as archive:
        for name in archive.namelist():
            if name.startswith('Contents/section') and name.endswith('.xml'):
                root=ET.fromstring(archive.read(name))
                for i,p in enumerate(root.iter(f'{{{HP}}}p')):
                    text=own_text(p)
                    if text:rows.append({'section':name,'index':i,'text':text})
    return rows

def _revise_package(source,output,replacements):
    hits=Counter()
    with zipfile.ZipFile(source) as original,zipfile.ZipFile(output,'w') as revised:
        for member in original.infolist():
            data=original.read(member.filename)
            if member.filename.startswith('Contents/section') and member.filename.endswith('.xml'):
                # Hancom requires its original namespace declarations even when
                # an XML serializer considers them unused. lxml preserves nsmap.
                root=ET.fromstring(data,ET.XMLParser(resolve_entities=False,no_network=True))
                changed=False
                for p in root.iter(f'{{{HP}}}p'):
                    text=own_text(p)
                    matching=[item for item in replacements if text==item['old']]
                    if len(matching)>1:raise ValueError('Ambiguous replacement specification')
                    if not matching:continue
                    item=matching[0];nodes=own_text_nodes(p)
                    if not nodes:raise ValueError('Paragraph has no own text nodes')
                    for n in nodes:
                        n.text=''
                        for child in list(n):n.remove(child)
                    nodes[0].text=item['new'];hits[item['id']]+=1;changed=True
                    # Cached source line positions can exceed a shortened new
                    # paragraph, making Hancom reject an otherwise valid package.
                    for cache in p.findall(f'./{{{HP}}}linesegarray'):
                        p.remove(cache)
                if changed:data=ET.tostring(root,encoding='UTF-8',xml_declaration=True,standalone=True)
            revised.writestr(member,data)
    for item in replacements:
        if hits[item['id']]!=item.get('expected',1):
            raise ValueError(f"{item['id']}: expected {item.get('expected',1)}, observed {hits[item['id']]}")
    return dict(hits)

def revise(source,output,replacements):
    source,output=Path(source),Path(output)
    if source.resolve()==output.resolve():
        raise ValueError('Source HWPX must remain unchanged; choose a separate output')
    descriptor,temporary=tempfile.mkstemp(prefix='.hwpx-review-',suffix='.tmp',dir=output.parent)
    os.close(descriptor)
    temporary=Path(temporary)
    try:
        hits=_revise_package(source,temporary,replacements)
        os.replace(temporary,output)
        return hits
    finally:
        temporary.unlink(missing_ok=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('input',type=Path);p.add_argument('--extract',type=Path);p.add_argument('--replacements',type=Path);p.add_argument('--output',type=Path)
    args=p.parse_args()
    if args.extract:
        rows=paragraphs(args.input);args.extract.write_text(json.dumps(rows,ensure_ascii=False,indent=2),encoding='utf-8')
        print(f"{len(rows)} text paragraphs; {sum('강근수 확인' in r['text'] for r in rows)} author markers")
    else:
        print(revise(args.input,args.output,json.loads(args.replacements.read_text(encoding='utf-8'))))
