"""Import an immutable heuristic literature report without promoting its claims.

The provider's exact-text checks are receipts, not a new semantic validation.
The separate snapshot database never writes to CTD or the running service.
"""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import sqlite3
from datetime import datetime, timezone
from contextlib import contextmanager


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


@contextmanager
def read_database(path):
    con = sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro', uri=True)
    con.row_factory = sqlite3.Row
    try:
        yield con
    finally:
        con.close()


def import_snapshot(report, database, *, expected_sha256, previous_manifest=None, deployment=None):
    report, database = Path(report), Path(database)
    report_bytes = report.read_bytes()
    if hashlib.sha256(report_bytes).hexdigest() != expected_sha256:
        raise ValueError('Frozen literature report SHA-256 mismatch')
    payload = json.loads(report_bytes)
    if payload.get('schema_version') != 'literature-report-v1':
        raise ValueError('Unsupported literature report schema')
    if database.exists():
        raise FileExistsError('Snapshot exists; choose a new versioned output')
    papers = payload['papers']
    ids = [str(p['pmid']) for p in papers]
    if any(not x.isdigit() for x in ids) or len(ids) != len(set(ids)):
        raise ValueError('Invalid or duplicate PMID')
    counts = payload['counts']
    if counts['papers_in_report'] != len(papers):
        raise ValueError('Report paper count mismatch')
    raw_assertions = sum(len(p.get('assertions', [])) for p in papers)
    if raw_assertions != counts['relation_candidates']:
        raise ValueError('Report assertion count mismatch')
    previous = set()
    if previous_manifest:
        old = json.loads(Path(previous_manifest).read_text(encoding='utf-8'))
        previous = {str(x) for topic in old['searches'].values() for x in topic['selected_pmids']}
        if len(previous) != old['unique_documents']:
            raise ValueError('Previous collection PMID count mismatch')
    database.parent.mkdir(parents=True, exist_ok=True)
    temporary = database.with_name(database.name+'.building')
    if temporary.exists():
        raise FileExistsError('Incomplete output already exists')
    con = sqlite3.connect(temporary)
    con.execute('PRAGMA foreign_keys=ON')
    con.executescript('''
      CREATE TABLE metadata(key TEXT PRIMARY KEY,value_json TEXT NOT NULL);
      CREATE TABLE paper(pmid TEXT PRIMARY KEY,payload_json TEXT NOT NULL);
      CREATE TABLE assertion(id TEXT PRIMARY KEY,pmid TEXT NOT NULL REFERENCES paper(pmid),
         subject TEXT NOT NULL,predicate TEXT NOT NULL,object TEXT NOT NULL,status TEXT NOT NULL,
         tier TEXT NOT NULL CHECK(tier='extracted'),payload_json TEXT NOT NULL);
      CREATE INDEX assertion_pmid ON assertion(pmid);
      CREATE INDEX assertion_status ON assertion(status);
    ''')
    imported = approved = stale = 0
    rules, dictionaries = set(), set()
    try:
        with con:
            for paper in papers:
                pmid = str(paper['pmid'])
                extraction = paper.get('extraction')
                current = bool(paper.get('extraction_current')) and not paper.get('retracted')
                if extraction:
                    if extraction['source_hash'] != extraction['source_text_hash']:
                        raise ValueError('Provider source hash mismatch: '+pmid)
                    rules.add(extraction['result']['rules_version'])
                    dictionaries.add(extraction['dictionary_version'])
                # Do not publish abstracts/quotations in Git; they stay in this local DB.
                entry = {k: paper[k] for k in ('pmid','doi','pmcid','title','pub_date','organism','oa_status','license','retracted','topics')}
                entry.update({'summary_current':current,'summary':extraction['result'].get('summary', {}) if extraction and current else None,
                    'context':extraction['result'].get('context', {}) if extraction and current else None,
                    'provider_provenance': {k:extraction[k] for k in ('extraction_id','source_id','source_hash','source_kind','source_url','source_license','engine_version','dictionary_version')} if extraction else None,
                    'source_report_sha256':expected_sha256,'also_in_original_1198':pmid in previous,
                    'interpretation':'heuristic summary; independent expert accuracy pending'})
                con.execute('INSERT INTO paper VALUES(?,?)',(pmid,canonical(entry)))
                for assertion in paper.get('assertions', []):
                    if not current or not assertion.get('source_span_valid'):
                        stale += 1
                        continue
                    if not 0 <= assertion['start'] < assertion['end'] or assertion['end']-assertion['start'] != len(assertion['quote']):
                        raise ValueError('Provider quotation/offset length mismatch')
                    status = assertion['display_status']
                    relation_payload = json.loads(assertion['payload'])
                    reviews = assertion.get('reviews', [])
                    original = {key:assertion[key] for key in ('subject','predicate','object')}
                    reviewed_relation = original
                    if status == 'approved':
                        latest = reviews[-1] if reviews else {}
                        if (latest.get('decision') not in ('approved','corrected') or
                            latest.get('fingerprint') != assertion['fingerprint'] or
                            not relation_payload.get('explicit_relation')):
                            raise ValueError('Approval without matching explicit review')
                        if latest['decision'] == 'corrected':
                            correction = latest.get('corrected_payload')
                            if isinstance(correction,str):correction=json.loads(correction)
                            if (not isinstance(correction,dict) or set(correction)!=set(original) or
                                any(not isinstance(x,str) or not x.strip() for x in correction.values())):
                                raise ValueError('Invalid reviewed correction')
                            reviewed_relation = correction
                        approved += 1
                    effective = assertion['effective']
                    if effective != reviewed_relation:
                        raise ValueError('Effective relation is not bound to its matching review')
                    item = dict(assertion)
                    item['payload'] = relation_payload
                    item.update({'pmid':pmid,'tier':'extracted','relation_is_causal':False,
                        'binding_label':False,'source_report_sha256':expected_sha256,
                        'provider_provenance':entry['provider_provenance'],
                        'span_validation':'provider frozen-report exact-text receipt; full source not revalidated by this importer'})
                    con.execute('INSERT INTO assertion VALUES(?,?,?,?,?,?,?,?)',
                        (assertion['assertion_id'],pmid,effective['subject'],effective['predicate'],effective['object'],status,'extracted',canonical(item)))
                    imported += 1
            metadata = {'imported_at':datetime.now(timezone.utc).isoformat(),
                'importer_source_sha256':digest(Path(__file__)),
                'provider_report_generated_at':payload['generated_at'],'source_report_sha256':expected_sha256,
                'source_report_path':str(report),'source_report_counts':counts,'papers':len(papers),
                'imported_current_assertions':imported,'approved_assertions':approved,'skipped_stale_assertions':stale,
                'rules_versions':sorted(rules),'dictionary_versions':sorted(dictionaries),
                'scientific_accuracy':payload['quality'],'source_namespaces_automatically_merged':False,
                'overlap':{'original_collection':len(previous),'prototype_collection':len(ids),
                    'shared_pmids':len(previous & set(ids)),'prototype_only_pmids':len(set(ids)-previous),
                    'union_unique_pmids':len(previous | set(ids))},
                'deployment_receipt_sha256':digest(deployment) if deployment else None,
                'prototype_code_reused':'source-preserved integrations/literature_service; consumes its versioned report schema',
                'same_pmid_not_independent_evidence':True,'CTD_graph_modified':False,
                'source_period_note':'Service selection follows its dated query; heterogeneous PubMed record dates are preserved, including updated book chapters'}
            if approved != counts['approved_relations']:
                raise ValueError('Approved relation count mismatch')
            con.execute('INSERT INTO metadata VALUES(?,?)',('snapshot',canonical(metadata)))
        if con.execute('PRAGMA integrity_check').fetchone()[0] != 'ok' or con.execute('PRAGMA foreign_key_check').fetchall():
            raise ValueError('Snapshot integrity failure')
        con.close()
        os.replace(temporary, database)
        metadata['snapshot_database_sha256'] = digest(database)
        return metadata
    except BaseException:
        con.close()
        temporary.unlink(missing_ok=True)
        raise


def status(database):
    with read_database(database) as con:
        return json.loads(con.execute("SELECT value_json FROM metadata WHERE key='snapshot'").fetchone()[0])


def search(database, *, pmid=None, query='', reviewed=None, limit=25, offset=0, assertions=False):
    if not 1 <= limit <= 500 or offset < 0:
        raise ValueError('Invalid pagination')
    clauses, args = [], []
    if pmid:
        clauses.append('pmid=?'); args.append(pmid)
    if reviewed is not None and assertions:
        clauses.append("status='approved'" if reviewed else "status!='approved'")
    if query:
        # Literal LIKE matching; wildcard input is not a different query language.
        escaped = query.replace('\\','\\\\').replace('%','\\%').replace('_','\\_')
        clauses.append("payload_json LIKE ? ESCAPE '\\'"); args.append('%'+escaped+'%')
    table, key = ('assertion','id') if assertions else ('paper','pmid')
    where = ' WHERE '+' AND '.join(clauses) if clauses else ''
    with read_database(database) as con:
        count = con.execute('SELECT count(*) FROM '+table+where,args).fetchone()[0]
        rows = con.execute('SELECT payload_json FROM '+table+where+' ORDER BY '+key+' LIMIT ? OFFSET ?', [*args,limit,offset]).fetchall()
    return {'items':[json.loads(r[0]) for r in rows],'total':count,'limit':limit,'offset':offset,
            'source_report_sha256':status(database)['source_report_sha256'],'research_only':True}
