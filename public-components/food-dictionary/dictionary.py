"""Versioned dictionary overlay. Source records and production databases are immutable."""
import argparse
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import unicodedata
import uuid

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
VERSION = 'food-concept-dictionary-v1.0.0'
RELEASE = 'food-concept-dictionary-20261005-0530'
NAMESPACE = uuid.UUID('738cf5b1-845e-4daf-b99b-4513dc848ea3')
PUBLIC = ROOT / 'fooddb_public_bulk_20261005_0350/public'
PRIVATE = ROOT / 'fooddb_enrichment_20261005_0134'
DELTA = ROOT / 'fooddb_enrichment_scope_20261005_0240'


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def identifier(kind, namespace, key):
    # Deliberately independent of labels, row order and release dates.
    return 'fcd:' + kind + ':' + str(uuid.uuid5(NAMESPACE, canonical([kind, namespace, key])))


def normalize(text):
    # Retrieval only: preserve stereo prefixes, punctuation, locants and preparation.
    return re.sub(r'\s+', ' ', unicodedata.normalize('NFC', text).strip()).casefold()


def retrieve(con, query, language='en', semantic_type=None, limit=10):
    """Candidate retrieval only. Homonyms and ambiguous forms retain their IDs."""
    if language not in ('ko','en') or not 1 <= limit <= 100:
        raise ValueError('Unsupported language or bounded page size')
    q=normalize(query)
    if not q:
        return []
    sql='SELECT l.concept_id,l.original_text,l.normalized_text,l.terminology_status,l.review_status,l.visibility,c.semantic_type FROM label l JOIN concept c USING(concept_id) WHERE l.language=?'
    args=[language]
    if semantic_type:
        sql+=' AND c.semantic_type=?';args.append(semantic_type)
    candidates={}
    for row in con.execute(sql+' ORDER BY l.label_id',args):
        cid,text,norm,status,review,visibility,kind=row
        score=1.0 if q==norm else SequenceMatcher(None,q,norm).ratio()
        if score<0.65:
            continue
        item={'concept_id':cid,'matched_original_label':text,'semantic_type':kind,'terminology_status':status,
              'review_status':review,'visibility':visibility,'retrieval_score':round(score,6),
              'mapping_status':'candidate','score_semantics':'Uncalibrated lexical retrieval; no equivalence or probability claim.'}
        if cid not in candidates or score>candidates[cid]['retrieval_score']:
            candidates[cid]=item
    return sorted(candidates.values(),key=lambda r:(-r['retrieval_score'],r['concept_id']))[:limit]


def decide_equivalence(left, right, evidence):
    if not evidence.get('source_verified') or not evidence.get('scope_review_passed'):
        return 'candidate', 'Evidence and scope verification required; names cannot prove equivalence.'
    if left['semantic_type'] != right['semantic_type']:
        return 'conflict', 'Different semantic types.'
    if left['semantic_type'] == 'molecule':
        for field in ('full_inchikey', 'original_inchi'):
            if not left.get(field) or left.get(field) != right.get(field):
                return 'conflict', 'Different or missing complete molecular identity.'
        if any(not x.get('single_component') or x.get('stereo_unspecified') for x in (left, right)):
            return 'candidate', 'Mixture or unspecified stereo requires a form decision.'
        if left.get('form_scope') != right.get('form_scope') or not left.get('form_scope'):
            return 'conflict', 'Salt, hydrate, tautomer or molecular form scope differs or is missing.'
        return 'source_verified', 'Complete structure and reviewed exact form agree.'
    if left['semantic_type'] == 'food_material':
        fields = ('authority_namespace', 'authority_id', 'organism', 'part', 'preparation', 'state')
        if not left.get('authority_namespace') or not left.get('authority_id'):
            return 'candidate', 'Food identity needs authority and explicit context.'
        if any(left.get(k) != right.get(k) for k in fields):
            return 'conflict', 'Native food identity or context differs.'
        if any(left.get(k) is None for k in fields[2:]):
            return 'candidate', 'Unresolved food facets cannot support a cross-source equivalence.'
        return 'source_verified', 'Native identity and all reviewed food facets agree.'
    return 'candidate', 'This semantic type requires a dedicated evidence rule.'


class Store:
    def __init__(self, path, inputs):
        if path.exists():
            raise FileExistsError('New immutable release required: ' + str(path))
        self.con = sqlite3.connect(path)
        self.con.executescript((HERE / 'schema.sql').read_text())
        self.con.execute('INSERT INTO release VALUES(?,?,?,?)', (RELEASE, VERSION, VERSION, canonical(inputs)))

    def evidence(self, key, uri, version, locator, payload, rights='private_pending', verified='source_verified'):
        eid = identifier('evidence', 'evidence', key)
        row = (eid, uri, version, locator, verified, rights, canonical(payload))
        existing = self.con.execute('SELECT * FROM evidence WHERE evidence_id=?', (eid,)).fetchone()
        if existing and existing != row:
            raise ValueError('Conflicting evidence identity')
        self.con.execute('INSERT OR IGNORE INTO evidence VALUES(?,?,?,?,?,?,?)', row)
        return eid

    def concept(self, namespace, key, kind, definition, facets, visibility='private'):
        cid = identifier('concept', namespace, key)
        row=(cid, namespace, key, kind, definition, 'project_operational_definition',
                          'Only the native record and explicitly stated distinguishing facets.',
                          'Other organisms, parts, processing states, chemical forms or measurement definitions.',
                          'Identity is source-qualified; a shared label alone establishes no equivalence.',
                          canonical(facets), visibility)
        prior=self.con.execute('SELECT * FROM concept WHERE concept_id=?',(cid,)).fetchone()
        if prior and prior!=row:
            raise ValueError('Conflicting concept identity or facets; explicit versioned review required')
        self.con.execute('INSERT OR IGNORE INTO concept VALUES(?,?,?,?,?,?,?,?,?,?,?)',row)
        return cid

    def label(self, cid, lang, text, status, evidence, visibility='private', role='alternate', review='source_verified', translator=None):
        if not text or not text.strip():
            return None
        lid = identifier('label', cid, canonical([lang, text, evidence, status]))
        preferred = self.con.execute("SELECT label_id,terminology_status FROM label WHERE concept_id=? AND language=? AND role='preferred'", (cid, lang)).fetchone()
        if role == 'preferred' and preferred and preferred[0] != lid:
            priority = {'official_source': 0, 'source_reported': 1, 'project_source_verified': 2, 'project_provisional': 3}
            if priority[status] < priority[preferred[1]]:
                self.con.execute("UPDATE label SET role='alternate' WHERE label_id=?", (preferred[0],))
            else:
                role = 'alternate'
        rights = 'project_authored_label_only' if status.startswith('project_') and visibility == 'public' else ('inherited_public_source' if visibility == 'public' else 'private_field_review_pending')
        self.con.execute('INSERT OR IGNORE INTO label VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)',
                         (lid, cid, lang, text, normalize(text), role, status, review, evidence, rights, visibility, None, translator))
        return lid

    def record(self, cid, namespace, native, version, raw, evidence, visibility):
        raw_json = canonical(raw)
        rid = identifier('record', namespace, canonical([native, version]))
        row = (rid, cid, namespace, str(native), version, raw_json, digest(raw_json.encode()), evidence, visibility)
        existing = self.con.execute('SELECT * FROM source_record WHERE record_id=?', (rid,)).fetchone()
        if existing and existing != row:
            raise ValueError('Conflicting duplicate native source record')
        self.con.execute('INSERT OR IGNORE INTO source_record VALUES(?,?,?,?,?,?,?,?,?)', row)
        return rid

    def review(self, kind, subject, obj, reason, payload):
        rid = identifier('review', kind, canonical([subject, obj, payload]))
        self.con.execute('INSERT OR IGNORE INTO review_item VALUES(?,?,?,?,?,?,?)',
                         (rid, kind, subject, obj, 'pending', reason, canonical(payload)))

    def mapping(self, subject, obj, predicate, status, evidence, decision, visibility='private', subject_version='MFDS2026-50', object_version='sealed-staged-projection-v1'):
        if predicate == 'skos:exactMatch' and (status != 'source_verified' or not decision.get('scope_review_passed') or not decision.get('identity_evidence_verified')):
            raise ValueError('Unverified exact match refused')
        if predicate=='skos:exactMatch':
            ev=self.con.execute('SELECT payload_json FROM evidence WHERE evidence_id=?',(evidence,)).fetchone()
            payload=json.loads(ev[0]) if ev else {}
            key=decision.get('full_inchikey')
            if not key or payload.get('full_inchikey')!=key or not payload.get('form_review_passed') or not payload.get('original_inchi_exact_match'):
                raise ValueError('Exact mapping requires pinned complete-key and original-form evidence')
            for cid in (subject,obj):
                row=self.con.execute('SELECT semantic_type,facets_json FROM concept WHERE concept_id=?',(cid,)).fetchone()
                if not row or row[0]!='molecule' or json.loads(row[1]).get('full_inchikey')!=key:
                    raise ValueError('Exact mapping endpoints must denote the same scoped molecule')
        mid = identifier('mapping', predicate, canonical([subject, obj, evidence]))
        values = (mid, subject, predicate, obj, status, decision['reason'], evidence, subject_version, object_version,
                  VERSION, 'Codex', None, '2026-10-05', visibility, canonical(decision))
        existing = self.con.execute('SELECT * FROM mapping WHERE mapping_id=?', (mid,)).fetchone()
        if existing and existing != values:
            raise ValueError('Mapping reimport conflicts with prior decision')
        self.con.execute('INSERT OR IGNORE INTO mapping VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', values)
        event = identifier('event', RELEASE, canonical(['add_mapping', mid]))
        self.con.execute('INSERT OR IGNORE INTO change_event VALUES(?,?,?,?,?,?,?)',
                         (event, RELEASE, 'add_mapping', mid, None, canonical({'status': status}), None))
        return mid

    def revoke(self, mapping_id, reason, reviewer):
        prior = self.con.execute('SELECT status,decision_json FROM mapping WHERE mapping_id=?', (mapping_id,)).fetchone()
        if not prior or not reviewer or not reason:
            raise ValueError('Existing mapping and explicit reviewer/reason required')
        self.con.execute("UPDATE mapping SET status='revoked' WHERE mapping_id=?", (mapping_id,))
        event = identifier('event', RELEASE, canonical(['revoke_mapping', mapping_id, reason, reviewer]))
        self.con.execute('INSERT OR IGNORE INTO change_event VALUES(?,?,?,?,?,?,?)',
                         (event, RELEASE, 'revoke_mapping', mapping_id, canonical(prior), canonical({'status': 'revoked', 'reason': reason}), reviewer))


def load_inputs():
    files = [PUBLIC / 'foods.ndjson.gz', PUBLIC / 'components.json', PUBLIC / 'compound_registry.json',
             PUBLIC / 'source_versions.json', PUBLIC / 'rights_manifest.json',
             PRIVATE / 'base_structures.json', PRIVATE / 'official_enrichment_seed.json', PRIVATE / 'translated_names.json',
             PRIVATE / 'cross_database_links.json', PRIVATE / 'food_bilingual_catalog.json',
             DELTA / 'delta_seed.private.json', DELTA / 'new_verified_names.private.json']
    result = {}
    for p in files:
        raw = p.read_bytes()
        result[str(p.relative_to(ROOT))] = {'sha256': digest(raw), 'bytes': len(raw)}
    for name in ('dictionary.py','schema.sql','evaluate.py'):
        raw=(HERE/name).read_bytes()
        result['pipeline/'+name]={'sha256':digest(raw),'bytes':len(raw)}
    # Assert sealed selected inputs against the earlier receipts, without hashing big databases.
    ledgers = [ROOT / ('deliverable_receipts_' + n + '.json') for n in ['20261005_0134','20261005_0240','20261005_0350']]
    pins = {}
    for ledger in ledgers:
        pins.update(json.loads(ledger.read_bytes())['files'])
    for name, rec in result.items():
        pin = pins.get(name)
        if pin and pin['sha256'] != rec['sha256']:
            raise ValueError('Sealed input changed: ' + name)
    with gzip.open(PUBLIC / 'foods.ndjson.gz', 'rt') as stream:
        foods = [json.loads(line) for line in stream]
    read = lambda p: json.loads(p.read_bytes())
    data = {'foods': foods, 'components': read(PUBLIC/'components.json'), 'registry': read(PUBLIC/'compound_registry.json'),
            'source_versions': read(PUBLIC/'source_versions.json'), 'base': read(PRIVATE/'base_structures.json'),
            'seed': read(PRIVATE/'official_enrichment_seed.json'), 'translations': read(PRIVATE/'translated_names.json'),
            'links': read(PRIVATE/'cross_database_links.json'), 'private_foods': read(PRIVATE/'food_bilingual_catalog.json'),
            'delta': read(DELTA/'delta_seed.private.json'), 'delta_names': read(DELTA/'new_verified_names.private.json')}
    # A sealed research projection, never an operational database; bounded alias extraction only.
    db = PRIVATE / 'enrichment_final.sqlite'
    con = sqlite3.connect(db.as_uri()+'?mode=ro&immutable=1', uri=True)
    rows = con.execute('SELECT full_inchikey,language,original_name FROM compound_alias ORDER BY full_inchikey,language,original_name LIMIT 5001').fetchall()
    con.close()
    if len(rows) > 5000:
        raise ValueError('Bounded alias extraction exceeded')
    data['aliases'] = rows
    result['bounded_sealed_alias_projection'] = {'rows': len(rows), 'sha256': digest(canonical(rows).encode()), 'origin': str(db), 'mode': 'readonly immutable; no active DB'}
    return data, result


def build(path):
    data, inputs = load_inputs()
    s = Store(path, inputs)
    compounds = {}
    sources = {r['source_id']: r for r in data['source_versions']['sources']}
    source_evidence = {key: s.evidence(key, r['source_url'], r['version'], r['archive_member'], r, 'CC0-1.0') for key, r in sources.items()}
    for food in data['foods']:
        facets = {k: food[k] for k in ('scientific_name_as_reported','part','preparation','basis_kind','source_id','fdc_id','data_type')}
        cid = s.concept('USDA:FDC:food', str(food['fdc_id']), 'food_material',
                        'A food material described as '+food['original_description']+' in native FDC record '+str(food['fdc_id'])+'; organism, part and processing remain only as reported.', facets, 'public')
        ev = source_evidence[food['source_id']]
        s.label(cid,'en',food['name_en'],'source_reported',ev,'public','preferred')
        s.record(cid,'USDA:FDC:food',food['fdc_id'],sources[food['source_id']]['version'],food,ev,'public')
    defs = {r['component_id']: r for r in data['components']['source_definitions']}
    ev = s.evidence('fdc-component-definitions', 'https://fdc.nal.usda.gov/download-datasets/', 'Native USDA definitions in three sealed releases', '/components.json/source_definitions', {'sha256':inputs[str((PUBLIC/'components.json').relative_to(ROOT))]['sha256']}, 'CC0-1.0')
    for r in data['components']['components']:
        d = defs[r['id']]
        kind = 'analytical_nutrient' if r['entity_type']=='nutrient' else 'analytical_measure'
        cid=s.concept('USDA:FDC:component',str(r['fdc_nutrient_id']),kind,
                      'An analytical reporting concept for '+d['name_en']+' in '+d['original_unit']+' under native USDA definition '+str(r['fdc_nutrient_id'])+'; molecular form remains unresolved unless separately supported.',
                      {**r,'original_unit':d['original_unit'],'native_definition':json.loads(d['source_definition_json'])},'public')
        s.label(cid,'en',d['name_en'],'source_reported',ev,'public','preferred')
        s.record(cid,'USDA:FDC:component',r['id'],'native-definition-v1',d,ev,'public')
    for mol in data['registry']['compounds']:
        key=mol['full_inchikey'];provenance=mol['source']
        ev=s.evidence(key,provenance['source_url'],'PubChem government-computed capture2026-10-05','Computed Descriptors',provenance,'NCBI-USGovernment-PublicDomain')
        cid=s.concept('InChI:complete',key,'molecule','A molecular entity with the complete structure identity '+key+', retaining source-reported stereo and form.',{'full_inchikey':key,'original_inchi':mol['inchi'],'pubchem_cid':mol['pubchem_cid']},'public')
        compounds[key]=cid
        s.label(cid,'en',mol['iupac_name'],'source_reported',ev,'public','preferred')
        s.record(cid,'PubChem:computed',mol['pubchem_cid'],'capture2026-10-05',mol,ev,'public')
    ev=s.evidence('staged-base-structures','local:sealed-research-projection','fooddb_enrichment_20261005_0134','base_structures.json',{'sha256':inputs[str((PRIVATE/'base_structures.json').relative_to(ROOT))]['sha256']})
    for r in data['base']:
        key=r['full_inchikey']
        if key not in compounds:
            compounds[key]=s.concept('InChI:complete',key,'molecule' if r['source_single_component'] else 'chemical_record_unresolved',
                                     'A staged chemical record distinguished by its complete structure representation; unresolved mixtures and stereo retain their uncertainty.',r)
        s.record(compounds[key],'staged:complete-structure',key,'projection2026-10-05',r,ev,'private')
    for key,lang,name in data['aliases']:
        if key in compounds and lang in ('en','ko'):
            s.label(compounds[key],lang,name,'source_reported',ev,role='preferred')
    term_concepts={}; identity_evidence={}
    for seed in (data['seed'],data['delta']):
        for e in seed['evidence']:
            identity_evidence[e['evidence_id']]=s.evidence(e['evidence_id'],e['source_uri'],e['source_version'],e.get('source_heading_id',e['evidence_id']),e)
        default=identity_evidence[seed['evidence'][0]['evidence_id']]
        for t in seed['terms']:
            # The delta adds evidence to an existing official term; never invent a second material.
            native=t.get('source_heading_id',t['term_id'])
            cid=s.concept('MFDS:official-heading',native,'regulatory_material_term',
                          'A regulatory material term distinguished by its official heading and stated specification scope; molecular identity is a separate assertion.',
                          {'cas_numbers':t['cas_numbers'],'formula_lines':t['formula_lines'],'paragraph_index':t['paragraph_index'],'section':t['section']})
            term_concepts[t['term_id']]=cid
            s.label(cid,'ko',t['original_ko'],'official_source',default,role='preferred')
            s.label(cid,'en',t['original_en'],'official_source',default,role='preferred')
            s.record(cid,'MFDS:official-heading',t['term_id'],t['source_version'],t,default,'private')
        terms={t['term_id']:t for t in seed['terms']}
        for b in seed['compound_bindings']:
            if b['review_state'] != 'verified' or b['match_grade'] != 'exact_form' or b['full_inchikey']!=b['matched_full_inchikey']:
                raise ValueError('Prior exact-form binding not verified')
            t=terms[b['term_id']]; key=b['full_inchikey']; cid=compounds[key]; ev=identity_evidence[b['identity_evidence_id']]
            endpoint=term_concepts[t['term_id']]
            facets=json.loads(s.con.execute('SELECT facets_json FROM concept WHERE concept_id=?',(endpoint,)).fetchone()[0])
            facets.update(full_inchikey=key,exact_form_scope_verified_by_prior_source_evidence=True)
            s.con.execute('UPDATE concept SET semantic_type=?,facets_json=?,definition=? WHERE concept_id=?',
                          ('molecule',canonical(facets),'A molecular entity denoted by the source heading '+t['original_en']+' with verified complete structure identity '+key+'; other forms are excluded.',endpoint))
            s.label(cid,'ko',t['original_ko'],'official_source',ev,role='preferred',review='prior_source_and_form_verified')
            s.label(cid,'en',t['original_en'],'official_source',ev,role='preferred',review='prior_source_and_form_verified')
            s.mapping(endpoint,cid,'skos:exactMatch','source_verified',ev,
                      {'scope_review_passed':True,'identity_evidence_verified':True,'full_inchikey':key,'reason':b['rationale'],'inherited_prior_machine_source_verification':True,'human_reviewed':False})
    for r in data['translations']:
        proof=json.loads(r['proof_json'])
        if not proof.get('form_review_passed') or not proof.get('independent_lexical_source_verified'):
            raise ValueError('Translation evidence is unverified')
        key=r['entity_id'];cid=compounds[key]
        ev=s.evidence(r['name_id'],r['source_uri'],r['source_version'],'Prior bilingual primary-study metadata search; exact form verification',proof,'project_authored_label_only')
        # Project-created short labels with citations only; no third-party abstracts, titles or synonyms in public output.
        visibility='public' if s.con.execute('SELECT visibility FROM concept WHERE concept_id=?',(cid,)).fetchone()[0]=='public' else 'private'
        s.label(cid,'ko',r['original_ko'],'project_source_verified',ev,visibility,'preferred','prior_primary_lexical_and_form_verified','Codex')
        s.label(cid,'en',r['original_en'],'project_source_verified',ev,visibility,review='prior_primary_lexical_and_form_verified',translator='Codex')
    food_concepts={}
    for r in data['private_foods']:
        ev=s.evidence('MFDS-food-'+r['id'],r['provenance']['source_url'],r['source_version'],r['id'],r['provenance'])
        cid=s.concept('MFDS:food-material-code',r['id'],'food_material','A source-described food material distinguished by organism, part and processing under its native material code.',{k:r.get(k) for k in ('taxon','part','preparation','state','authority_namespace','authority_id')})
        food_concepts[r['id']]=cid
        s.record(cid,'MFDS:food-material-code',r['id'],r['source_version'],r,ev,'private')
        for name in r['ko_names']:s.label(cid,'ko',name,'official_source',ev,role='preferred')
        for name in r['en_names']:
            generated=name==r.get('translated_descriptor')
            s.label(cid,'en',name,'project_provisional' if generated else 'official_source',ev,role='alternate',review='scope_conflict' if r.get('official_english_scope_conflict') else 'source_reported',translator='Codex' if generated else None)
        s.review('food_scope',cid,None,'Organism-level English names and translated food descriptors do not prove processed-food equivalence.',r)
    for r in data['links']:
        if r['status']=='identity_verified':continue
        subject=term_concepts.get(r['source_id']) if r['entity_type']=='compound' else food_concepts.get(r['source_id'])
        obj=compounds.get(r['target_id'])
        s.review('prior_mapping_'+r['status'],subject,obj,r['review_reason'],r)
    normalized=defaultdict(set)
    for cid,lang,text in s.con.execute('SELECT concept_id,language,normalized_text FROM label'):
        normalized[(lang,text)].add(cid)
    for (lang,text), ids in sorted(normalized.items()):
        if len(ids)>1:
            s.review('lexical_collision',None,None,'Shared normalized label is a retrieval candidate, never an automatic merge.',{'language':lang,'normalized_text':text,'concept_ids':sorted(ids)})
    s.con.commit()
    if s.con.execute('PRAGMA foreign_key_check').fetchall() or s.con.execute('PRAGMA integrity_check').fetchone()[0]!='ok':
        raise ValueError('Dictionary integrity failed')
    export(s.con, path.parent, inputs)
    s.con.close()
    return path


def export(con, folder, inputs):
    con.row_factory=sqlite3.Row
    public=folder/'public';public.mkdir()
    def rows(table, where='', args=()):
        return [dict(r) for r in con.execute('SELECT * FROM '+table+(' WHERE '+where if where else '')+' ORDER BY 1',args)]
    all_concepts=rows('concept'); all_labels=rows('label'); all_maps=rows('mapping'); all_reviews=rows('review_item')
    public_ids={r['concept_id'] for r in all_concepts if r['visibility']=='public'}
    labels=[r for r in all_labels if r['visibility']=='public' and r['concept_id'] in public_ids]
    # Private preferred terms do not hide the best public source label.
    grouped=defaultdict(list)
    for r in labels:grouped[(r['concept_id'],r['language'])].append(r)
    for group in grouped.values():
        ordered=sorted(group,key=lambda x:(x['role']!='preferred',x['terminology_status']=='project_provisional',x['label_id']))
        for i,r in enumerate(ordered):r['role']='preferred' if i==0 else 'alternate'
    preferences=defaultdict(lambda:{'ko':None,'en':None})
    for r in labels:
        if r['role']=='preferred':preferences[r['concept_id']][r['language']]=r['original_text']
    selected_evidence={r['evidence_id'] for r in labels}
    records=rows('source_record',"visibility='public'")
    selected_evidence.update(r['evidence_id'] for r in records)
    ev=[]
    for r in rows('evidence'):
        if r['evidence_id'] in selected_evidence:
            # Minimized public provenance, with the private evidence payload retained only in research DB.
            ev.append({k:r[k] for k in ('evidence_id','source_uri','source_version','locator','verification_status','rights_status')})
    datasets={'concepts.ndjson':[{**r,'preferred_labels':preferences[r['concept_id']]} for r in all_concepts if r['concept_id'] in public_ids],
              'labels.ndjson':labels,'source_crosswalk.ndjson':[{k:r[k] for k in ('record_id','concept_id','source_namespace','native_id','source_version','raw_record_sha256','evidence_id')} for r in records],
              'evidence.ndjson':ev,'mappings.ndjson':[r for r in all_maps if r['visibility']=='public' and r['subject_id'] in public_ids and r['object_id'] in public_ids]}
    for name, values in datasets.items():
        (public/name).write_text(''.join(canonical(v)+'\n' for v in values))
    def coverage(concepts, labels):
        counts=Counter(r['semantic_type'] for r in concepts); out={}
        for kind,n in sorted(counts.items()):
            ids={r['concept_id'] for r in concepts if r['semantic_type']==kind}
            allowed={'source_verified','prior_source_and_form_verified','prior_primary_lexical_and_form_verified'}
            ko={r['concept_id'] for r in labels if r['concept_id'] in ids and r['language']=='ko' and r['review_status'] in allowed and r['terminology_status']!='project_provisional'}
            en={r['concept_id'] for r in labels if r['concept_id'] in ids and r['language']=='en' and r['review_status'] in allowed and r['terminology_status']!='project_provisional'}
            out[kind]={'denominator':n,'source_verified_ko_labels':len(ko),'verified_bilingual_concepts':len(ko&en),'missing_ko_labels':n-len(ko)}
        return out
    result={'release_id':RELEASE,'schema_version':VERSION,'counts':{'concepts':len(all_concepts),'labels':len(all_labels),'source_records':con.execute('SELECT count(*) FROM source_record').fetchone()[0],'source_verified_exact_mappings':sum(r['status']=='source_verified' and r['predicate_id']=='skos:exactMatch' for r in all_maps),'review_items':len(all_reviews),'human_reviewed_labels':0},
            'public_counts':{name:len(values) for name,values in datasets.items()},'research_coverage':coverage(all_concepts,all_labels),'public_coverage':coverage(datasets['concepts.ndjson'],labels),
            'canonical_chemical_structure_denominator':len({r['identity_key'] for r in all_concepts if r['identity_namespace']=='InChI:complete'}),
            'source_records_are_not_unique_food_species':True,'name_only_equivalence_merges':0,'inferred_transitive_merges':0,'nutrient_to_molecule_mappings':0,'human_gold_precision_recall':None,'SOTA_superiority_evaluated':False}
    canonical_compounds=[r for r in all_concepts if r['identity_namespace']=='InChI:complete']
    result['canonical_structure_bilingual_coverage']=coverage(canonical_compounds,all_labels)
    source_molecules={r['concept_id'] for r in all_concepts if r['semantic_type']=='molecule' and r['identity_namespace']!='InChI:complete'}
    result['counts']['source_molecule_concepts_with_exact_link']=len({
        endpoint for r in all_maps
        if r['status']=='source_verified' and r['predicate_id']=='skos:exactMatch'
        for endpoint in (r['subject_id'],r['object_id']) if endpoint in source_molecules})
    (folder/'RESULT.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    (public/'coverage.json').write_text(json.dumps({k:result[k] for k in ('release_id','schema_version','public_counts','public_coverage','name_only_equivalence_merges','inferred_transitive_merges','nutrient_to_molecule_mappings','human_gold_precision_recall')},ensure_ascii=False,indent=2)+'\n')
    (folder/'conflict_review.private.ndjson').write_text(''.join(canonical(r)+'\n' for r in all_reviews))
    # Source-qualified, concept-level deterministic split: labels never leak across train/test concepts.
    strata=defaultdict(list)
    for r in sorted(all_reviews,key=lambda r:r['review_id']):strata[r['kind']].append(r)
    held=[]
    while len(held)<200 and any(strata.values()):
        for kind in sorted(strata):
            if strata[kind] and len(held)<200:held.append(strata[kind].pop(0))
    cases=[{'review_id':r['review_id'],'subject_id':r['subject_id'],'object_id':r['object_id'],'sampling_frame':r['kind'],'related_concept_ids':json.loads(r['payload_json']).get('concept_ids',[])} for r in held]
    cases.extend({'review_id':r['mapping_id'],'subject_id':r['subject_id'],'object_id':r['object_id'],'sampling_frame':'prior_source_verified_exact_mapping','related_concept_ids':[]} for r in all_maps)
    parents={}
    def find(x):
        parents.setdefault(x,x)
        while parents[x]!=x:
            parents[x]=parents[parents[x]];x=parents[x]
        return x
    for r in cases:
        ids=sorted(set([x for x in (r['subject_id'],r['object_id']) if x]+r['related_concept_ids']))
        for item in ids[1:]:parents[find(item)]=find(ids[0])
    evaluation=[]
    for r in cases:
        ids=[x for x in (r['subject_id'],r['object_id']) if x]+r['related_concept_ids']
        family=find(ids[0]) if ids else r['review_id']
        split='held_out' if int(digest(family.encode())[:8],16)%5==0 else 'development'
        evaluation.append({**r,'concept_group_id':family,'split':split,'reviewer':None,'gold_predicate':None,'gold_equivalent':None,'prediction_equivalent':None,'retrieved_object_ids':[],'status':'awaiting_actual_review'})
    (folder/'evaluation_review_set.private.ndjson').write_text(''.join(canonical(r)+'\n' for r in evaluation))
    rights={'release_id':RELEASE,'publication_eligible':True,'field_decisions':[
        {'fields':'USDA source labels and source-qualified food/component crosswalks','license':'CC0-1.0','source':'https://fdc.nal.usda.gov/'},
        {'fields':'PubChem53 computed structure identities and IUPAC labels','license':'LicenseRef-NCBI-USGovernment-PublicDomain','source':'https://www.ncbi.nlm.nih.gov/home/about/policies/','restriction':'Government-computed capture only; no contributed synonym catalog.'},
        {'fields':'Five short project-created Korean transliterations and corresponding English aliases','license':'LicenseRef-ProjectGeneratedLabels','basis':'Project-authored terminology with primary lexical citations; no third-party title, abstract, protocol or description reproduced. Source verification is lexical/form verification, not a human review.'},
        {'fields':'Internal concept IDs and project operational definitions','license':'LicenseRef-ProjectGeneratedDictionary','basis':'New project-generated material in authorized companion export.'}],
        'excluded':['All MFDS headings and aliases pending field-specific reuse review','All private FooDB aliases and staged structure payloads','All private source records, mappings, reviews and evaluation labels','RDA file/API data not acquired','FoodOn/ChEBI/LanguaL terms not imported'],
        'attribution':'USDA FoodData Central and NCBI PubChem; individual primary lexical-source URLs are in evidence.ndjson. No source agency endorsement.',
        'prior_public_releases_unchanged':True,'site_activation_executed':False}
    (public/'rights_manifest.json').write_text(json.dumps(rights,ensure_ascii=False,indent=2)+'\n')
    (public/'contract.json').write_text(json.dumps({'schema_version':VERSION,'identity_rule':'Stable UUIDv5 over namespace and exact native/complete structure identity, independent of label and release.','normalization':'NFC, whitespace folding and casefold for retrieval only; original labels preserved.','source_crosswalk':'Reversible native identifiers and versioned source-record hashes.','relations':'SKOS exact/close/broad/narrow/related; exact requires verified exact-form scope. No transitive union or identity rewriting.','definition_status':'Project operational definitions; no human curator acceptance inferred.','unknown_ko':'Null/absent Korean labels must remain visible as untranslated.','integration':'Immutable companion only. Site owner must adapt contract and approve scoped activation; do not replace production IDs or mutate source records.'},ensure_ascii=False,indent=2)+'\n')
    manifest={'release_id':RELEASE,'schema_version':VERSION,'publication_eligible':True,'files':{p.name:{'bytes':p.stat().st_size,'sha256':digest(p.read_bytes())} for p in sorted(public.iterdir())},'selected_input_pins':{k:v for k,v in inputs.items() if k.startswith('fooddb_public_bulk_')},'pipeline_pins':{k:v for k,v in inputs.items() if k.startswith('pipeline/')},'counts':result['public_counts']}
    (public/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE/'research_dictionary.sqlite')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    build(args.output)
    print((args.output.parent/'RESULT.json').read_text())
