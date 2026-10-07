"""Native ChemProt gold, source offsets and a separately sealed evaluation label file."""
from __future__ import annotations

from collections import Counter
import hashlib
import io
import json
from pathlib import Path
import random
import zipfile

from .extraction import load_drugprot

CPR_LABELS = ("CPR:3", "CPR:4", "CPR:5", "CPR:6", "CPR:9")


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def mark_context(text, chemical, gene, flank=240):
    """Source-derived pair context. Entity boundaries and types, never relation labels."""
    if not all(0 <= e['start'] < e['end'] <= len(text) and text[e['start']:e['end']] == e['text'] for e in (chemical, gene)):
        raise ValueError('Entity source offset mismatch')
    a, b = sorted((chemical, gene), key=lambda e:e['start'])
    if a['end'] > b['start']:
        raise ValueError('Overlapping chemical/gene mentions cannot be represented by this protocol')
    start, end = max(0, a['start']-flank), min(len(text), b['end']+flank)
    if b['start']-a['end'] > 1200:
        chunks = [(max(0,e['start']-flank),min(len(text),e['end']+flank)) for e in (a,b)]
    else:
        chunks = [(start,end)]
    rendered = []
    for lo,hi in chunks:
        fragment = text[lo:hi]
        for e, tag in sorted(((chemical,'CHEM'),(gene,'GENE')),key=lambda x:x[0]['start'],reverse=True):
            if lo <= e['start'] and e['end'] <= hi:
                i,j = e['start']-lo,e['end']-lo
                fragment = fragment[:i] + f'[{tag}] ' + fragment[i:j] + f' [/{tag}]' + fragment[j:]
        rendered.append(fragment)
    return ' [SEP] '.join(rendered)


def pair_rows(documents):
    for d in documents:
        for a in d['entities']:
            if a['type'] != 'CHEMICAL':
                continue
            for b in d['entities']:
                if not b['type'].startswith('GENE'):
                    continue
                yield {'pmid':d['pmid'],'chemical_id':a['id'],'gene_id':b['id'],
                       'text':mark_context(d['text'],a,b)}


def read_cpr_gold(path):
    gold = {}
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        pmid,label,a,b = line.split('\t')
        if label not in CPR_LABELS or not a.startswith('Arg1:') or not b.startswith('Arg2:'):
            raise ValueError('Unsupported native ChemProt gold row')
        gold.setdefault((pmid,a[5:],b[5:]),set()).add(label)
    return gold


def prepare_chemprot(archive, output, exposed_ids=(), exposed_text_hashes=(), seed=20261003):
    """Only training/development gold are parsed. Test labels are copied and hashed."""
    out = Path(output)
    out.mkdir(parents=True,exist_ok=False)
    exposed = set(map(str,exposed_ids))
    seen_hashes = set(exposed_text_hashes)
    ledger, documents_by_split = {}, {}
    with zipfile.ZipFile(archive) as outer:
        for split,native in [('train','training'),('dev','development'),('test','test_gs')]:
            name = f'ChemProt_Corpus/chemprot_{native}.zip'
            with zipfile.ZipFile(io.BytesIO(outer.read(name))) as inner:
                files = {}
                for suffix in ['abstracts','entities','gold_standard']:
                    matches = [n for n in inner.namelist() if '/._' not in n and n.endswith('.tsv') and f'_{suffix}' in Path(n).name]
                    if len(matches) != 1:
                        raise ValueError('Ambiguous ChemProt member')
                    target = out/f'{split}_{suffix}.tsv'
                    target.write_bytes(inner.read(matches[0]))
                    files[suffix] = target
                documents = load_drugprot(files['abstracts'],files['entities'])
                excluded = [d['pmid'] for d in documents if d['pmid'] in exposed or hashlib.sha256(d['text'].encode()).hexdigest() in seen_hashes]
                retained = [d for d in documents if d['pmid'] not in set(excluded)]
                documents_by_split[split] = retained
                rows = list(pair_rows(retained))
                if split != 'test':
                    gold = read_cpr_gold(files['gold_standard'])
                    for row in rows:
                        labels = gold.get((row['pmid'],row['chemical_id'],row['gene_id']),set())
                        row['labels'] = [int(label in labels) for label in CPR_LABELS]
                    if split == 'train':
                        positives = [r for r in rows if any(r['labels'])]
                        negatives = [r for r in rows if not any(r['labels'])]
                        random.Random(seed).shuffle(negatives)
                        rows = positives + negatives[:3*len(positives)]
                        random.Random(seed).shuffle(rows)
                with (out/f'{split}.jsonl').open('w',encoding='utf-8') as f:
                    for row in rows:
                        f.write(json.dumps(row,ensure_ascii=False)+'\n')
                ledger[split] = {'original_documents':len(documents),'retained_documents':len(retained),
                                 'excluded_previously_exposed_pmids':sorted(excluded),'candidate_pairs':len(rows),
                                 'pmids':[d['pmid'] for d in retained],
                                 'text_hashes':[hashlib.sha256(d['text'].encode()).hexdigest() for d in retained],
                                 'gold_sha256':sha(files['gold_standard']),
                                 'inputs_sha256':sha(out/f'{split}.jsonl')}
    for a,b in [('train','dev'),('train','test'),('dev','test')]:
        if set(ledger[a]['pmids']) & set(ledger[b]['pmids']) or set(ledger[a]['text_hashes']) & set(ledger[b]['text_hashes']):
            raise ValueError('Official ChemProt evaluation boundary contains source overlap')
    protocol = {'protocol_id':'chemprot-biomedbert-20261003-v1','seed':seed,'labels':CPR_LABELS,
                'endpoint':'Presence of a native shared-task CPR:3/4/5/6/9 relation for a supplied chemical–gene mention pair',
                'score':'1-product(1-sigmoid(label_logits))','gold_mentions_used':True,
                'clinical_or_food_efficacy_evaluated':False,'test_gold_parsed_during_preparation':False,
                'negative_sampling':'train only, at most 3 negatives per positive; dev/test retain all mention pairs',
                'exposure_audit':'Previously used DrugProt/BioRED and public-corpus PMIDs/text hashes excluded from every partition',
                'raw_bert_pretraining_overlap':'unknown; PubMed/PMC pretraining is not annotation-level independence',
                'archive_sha256':sha(archive),'splits':ledger}
    write_json(out/'protocol.json',protocol)
    return protocol


def prepare_bioredirect(archive, output, parser, exposed_ids=(), exposed_text_hashes=(), seed=20261003):
    """Train/dev are old official partitions; BC8 400-paper gold stays sealed."""
    out = Path(output)
    out.mkdir(parents=True,exist_ok=False)
    exposed, seen = set(map(str,exposed_ids)),set(exposed_text_hashes)
    labels = ['Association','Bind','Comparison','Conversion','Cotreatment','Drug_Interaction',
              'Negative_Correlation','Positive_Correlation']
    splits = {}
    with zipfile.ZipFile(archive) as z:
        for split,member in [('train','bioredirect_train.pubtator'),('dev','bioredirect_dev.pubtator'),
                             ('test','bioredirect_bc8_test.pubtator')]:
            payload = z.read('bioredirect/'+member)
            # Remove native relation lines BEFORE invoking the parser for test inputs.
            text = payload.decode()
            stripped = '\n'.join(line for line in text.splitlines() if '|' in line.split('\t')[0] or not line.strip()
                                 or (len(line.split('\t'))==6 and line.split('\t')[1].isdigit()))
            documents = parser(stripped)
            excluded = []
            if split=='test':
                excluded = [d['pmid'] for d in documents if d['pmid'] in exposed or hashlib.sha256(d['text'].encode()).hexdigest() in seen]
                documents = [d for d in documents if d['pmid'] not in set(excluded)]
            lookup = {}
            if split!='test':
                lookup = {d['pmid']:d for d in parser(text)}
            rows, native_entities = [],[]
            for d in documents:
                groups = {}
                for e in d['entities']:
                    if e['kind'] not in {'gene','component'}:
                        continue
                    for identifier in e['native_ids']:
                        groups.setdefault((e['kind'],identifier),[]).append(e)
                # One canonical entity pair per paper, with closest non-overlapping mentions.
                for (kind,a),chemicals in groups.items():
                    if kind!='component':continue
                    for (other,b),genes in groups.items():
                        if other!='gene':continue
                        options = [(c,g) for c in chemicals for g in genes if c['end']<=g['start'] or g['end']<=c['start']]
                        if not options:raise ValueError('No non-overlapping canonical mentions')
                        c,g = min(options,key=lambda pair:abs(pair[0]['start']-pair[1]['start']))
                        row = {'pmid':d['pmid'],'chemical_id':a,'gene_id':b,'text':mark_context(d['text'],c,g)}
                        if split!='test':
                            positive = {r['predicate'] for r in lookup[d['pmid']]['relations'] if {r['entity1'],r['entity2']}=={a,b}}
                            if positive-set(labels):raise ValueError('Unrecognized native relation type')
                            row['labels'] = [int(l in positive) for l in labels]
                        rows.append(row)
                native_entities.append(d)
            # Raw test gold bytes are retained separately; score parses them after model lock.
            (out/f'{split}_native.pubtator').write_bytes(payload)
            with (out/f'{split}.jsonl').open('w',encoding='utf-8') as f:
                for row in rows:f.write(json.dumps(row,ensure_ascii=False)+'\n')
            write_json(out/f'{split}_documents.json',native_entities)
            splits[split]={'retained_documents':len(documents),'excluded_previously_exposed_pmids':excluded,
                           'candidate_pairs':len(rows),'inputs_sha256':sha(out/f'{split}.jsonl'),
                           'gold_sha256':sha(out/f'{split}_native.pubtator'),'pmids':[d['pmid'] for d in documents],
                           'text_hashes':[hashlib.sha256(d['text'].encode()).hexdigest() for d in documents]}
    for a,b in [('train','dev'),('train','test'),('dev','test')]:
        if set(splits[a]['pmids']) & set(splits[b]['pmids']) or set(splits[a]['text_hashes']) & set(splits[b]['text_hashes']):
            raise ValueError('BioRED evaluation boundary overlap')
    protocol={'protocol_id':'biored-bc8-biomedbert-20261003-v1','seed':seed,'labels':labels,
              'endpoint':'A published native BioRED relation for a supplied canonical chemical–gene pair',
              'score':'1-product(1-sigmoid(label_logits))','gold_mentions_used':True,
              'clinical_or_food_efficacy_evaluated':False,'test_gold_parsed_during_preparation':False,
              'negative_sampling':'none; all canonical chemical–gene pairs per document retained',
              'raw_bert_pretraining_overlap':'unknown; PubMed/PMC pretraining does not establish text-level independence',
              'archive_sha256':sha(archive),'splits':splits}
    write_json(out/'protocol.json',protocol)
    return protocol
