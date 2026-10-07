"""Owned isolated food DB migration; never opens an operational source DB."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import unicodedata

HERE = Path(__file__).resolve().parent
STRICT = 'strict_positive_source_named_single_analyte'
FULL_KEY = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')

def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)

def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()

def language(name):
    name = unicodedata.normalize('NFC', name)
    if re.search('[가-힣]', name): return 'ko'
    if re.search('[A-Za-z]', name): return 'en'
    return 'und'

def connect(path):
    con = sqlite3.connect(path)
    con.execute('PRAGMA foreign_keys=ON')
    con.execute('PRAGMA busy_timeout=2000')
    return con

def initialize(con):
    sql = (HERE / 'schema.sql').read_text()
    con.executescript(sql)
    def normalized(value): return ' '.join(value.lower().replace('if not exists ', '').split()).rstrip(';')
    for kind, ending in [('TRIGGER', 'END;'), ('VIEW', ';')]:
        for definition in re.finditer(r'CREATE ' + kind + r' IF NOT EXISTS ([a-z_]+)\b.*?' + ending, sql, re.S):
            name, desired = definition.group(1), definition.group(0)
            existing = con.execute('SELECT sql FROM sqlite_master WHERE type=? AND name=?', (kind.lower(), name)).fetchone()[0]
            if normalized(existing) != normalized(desired):
                con.execute('DROP ' + kind + ' ' + name)
                con.executescript(desired)

def source_identity(inputs):
    return digest(sorted(({'sha256': r['sha256'], 'bytes': r['bytes'],
                           'table': 'standards' if r['remote_path'].endswith('/strict_verified_food_compounds.parquet') else 'audit'} for r in inputs), key=lambda r: r['table']))

def import_projection(con, path, expected_sha256):
    actual_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_sha != expected_sha256: raise ValueError('projection SHA mismatch')
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        header = json.loads(next(stream))
        if header['kind'] != 'header': raise ValueError('missing projection header')
        snapshot = source_identity(header['inputs'])
        existing = con.execute('SELECT import_complete,projection_sha256 FROM source_snapshot WHERE snapshot_id=?', (snapshot,)).fetchone()
        if existing and existing[0]:
            if existing[1] != actual_sha: raise ValueError('same input snapshot but different projection; require explicit reconciliation')
            return {'snapshot_id': snapshot, 'already_imported': True, 'new_raw_rows': 0}
        caches = {'compound': {}, 'representation': set(), 'identifier': set(), 'alias': set(), 'food': set(), 'food_identifier': set(), 'context': set()}
        counts = {'standards': 0, 'audit': 0}
        with con:
            con.execute('INSERT INTO source_snapshot VALUES(?,?,?,?,?,?,?,0)',
                        (snapshot, 'food_compound_extension_source_analyte_audit_v2', '2026-10-03T08:04:36.601855+00:00',
                         'https://github.com/SWi1/FooDB_polyphenol_analysis', '3d2bcf9523911fabc8a07c731226ac87a2af262e', canonical(header['inputs']), actual_sha))
            for item in header['inputs']:
                table = 'standards' if item['remote_path'].endswith('/strict_verified_food_compounds.parquet') else 'audit'
                con.execute('INSERT INTO source_file VALUES(?,?,?,?,?)', (snapshot, table, item['sha256'], item['bytes'], 0))
            footer = None
            for line in stream:
                item = json.loads(line)
                if item['kind'] == 'footer':
                    if footer is not None: raise ValueError('duplicate footer')
                    footer = item; continue
                if footer is not None: raise ValueError('row after footer')
                table, row = item['table'], item['record']
                if table not in counts: raise ValueError('unknown source table')
                counts[table] += 1
                raw = canonical(row)
                con.execute('INSERT INTO raw_record VALUES(?,?,?,?,?)', (snapshot, table, row['id'], hashlib.sha256(raw.encode()).hexdigest(), raw))
                key = row['full_inchikey']
                if not FULL_KEY.fullmatch(key): raise ValueError('invalid full InChIKey')
                flags=(int(row['standard_structure_verified']), int(row['single_component']), int(row['standard_unspecified_stereo']))
                prior_flags=caches['compound'].get(key)
                if prior_flags is None:
                    existing_flags=con.execute('SELECT source_structure_verified,source_single_component,source_stereo_unspecified FROM compound WHERE full_inchikey=?',(key,)).fetchone()
                    prior_flags=tuple(existing_flags) if existing_flags is not None else None
                if prior_flags is not None and prior_flags!=flags:
                    raise ValueError('Conflicting structure flags for a complete key; explicit source reconciliation required')
                if key not in caches['compound']:
                    con.execute('INSERT INTO compound VALUES(?,?,?,?) ON CONFLICT(full_inchikey) DO NOTHING',(key,*flags))
                    caches['compound'][key]=flags
                rep = (key, row.get('source_inchi'), row.get('canonical_isomeric_smiles'))
                if rep not in caches['representation']:
                    con.execute('INSERT INTO structure_representation VALUES(?,?,?,?,?)', (digest((snapshot, *rep)), key, snapshot, rep[1], rep[2])); caches['representation'].add(rep)
                for namespace, field in [('FooDB-public-compound', 'public_compound_id'), ('FooDB-native-compound', 'source_id')]:
                    value = row.get(field)
                    if value not in [None, '']:
                        token = (snapshot, namespace, str(value), key)
                        if token not in caches['identifier']:
                            con.execute('INSERT INTO compound_identifier_assertion VALUES(?,?,?,?)', token); caches['identifier'].add(token)
                alias_fields = ['name'] + (['orig_source_name'] if row['strict_status'] == STRICT else [])
                for field in alias_fields:
                    name = row.get(field)
                    if name not in [None, '']:
                        token = (key, name, language(name), snapshot)
                        if token not in caches['alias']:
                            con.execute('INSERT INTO compound_alias(full_inchikey,original_name,search_name,language,snapshot_id) VALUES(?,?,?,?,?)',
                                        (key, name, unicodedata.normalize('NFC', name).casefold(), token[2], snapshot)); caches['alias'].add(token)
                if table != 'audit': continue
                food_id = row.get('food_id')
                if food_id in [None, '']: raise ValueError('missing native food ID')
                if food_id not in caches['food']:
                    con.execute('INSERT INTO food VALUES(?,?)', (snapshot, food_id)); caches['food'].add(food_id)
                public = row.get('public_food_id')
                if public not in [None, '']:
                    token = (snapshot, food_id, 'FooDB-public-food', public)
                    if token not in caches['food_identifier']:
                        con.execute('INSERT INTO food_identifier_assertion VALUES(?,?,?,?)', token); caches['food_identifier'].add(token)
                context = (snapshot, food_id, row.get('orig_food_common_name'), row.get('orig_food_part'))
                context_id = digest(context)
                if context_id not in caches['context']:
                    con.execute('INSERT INTO food_context(context_id,snapshot_id,native_food_id,original_common_name,original_part) VALUES(?,?,?,?,?)', (context_id, *context)); caches['context'].add(context_id)
                state = 'accepted_strict' if row['strict_status'] == STRICT else 'retained_held'
                unit = 'mg_per_100g' if row.get('orig_unit') == 'mg/100g' else None
                con.execute('INSERT INTO observation(snapshot_id,source_row_id,full_inchikey,context_id,positive_value,original_unit,normalized_unit_code,original_method,original_analyte_id,original_analyte_name,source_citation,source_citation_type,source_pmid,strict_status,quality_state) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                            (snapshot, row['id'], key, context_id, row.get('positive_value'), row.get('orig_unit'), unit,
                             row.get('orig_method'), row.get('orig_source_id'), row.get('orig_source_name'), row.get('citation'), row.get('citation_type'), row.get('source_pmid'), row['strict_status'], state))
            if footer is None or footer['counts'] != counts or footer['total_rows'] != sum(counts.values()): raise ValueError('footer counts mismatch')
            for table, n in counts.items(): con.execute('UPDATE source_file SET row_count=? WHERE snapshot_id=? AND source_table=?', (n, snapshot, table))
            con.execute('UPDATE source_snapshot SET import_complete=1 WHERE snapshot_id=?', (snapshot,))
            if con.execute('PRAGMA foreign_key_check').fetchall(): raise ValueError('foreign-key violation')
        return {'snapshot_id': snapshot, 'already_imported': False, 'new_raw_rows': sum(counts.values()), 'source_counts': counts}

def add_term_seed(con, seed):
    with con:
        for evidence in seed['evidence']:
            body = canonical(evidence); token = digest(evidence)
            existing = con.execute('SELECT evidence_sha256 FROM authority_evidence WHERE evidence_id=?', (evidence['evidence_id'],)).fetchone()
            if existing and existing[0] != token: raise ValueError('existing authority evidence differs; register a new evidence version')
            con.execute('INSERT OR IGNORE INTO authority_evidence VALUES(?,?,?,?,?,?,?,?)',
                        (evidence['evidence_id'], evidence['authority'], evidence['source_uri'], evidence['source_version'], evidence['retrieved_utc'], evidence['verification_method'], token, body))
        for term in seed['terms']:
            existing = con.execute('SELECT raw_source_json FROM official_term WHERE term_id=?', (term['term_id'],)).fetchone()
            if existing and existing[0] != canonical(term): raise ValueError('existing official term differs; register a new term version')
            con.execute('INSERT OR IGNORE INTO official_term VALUES(?,?,?,?,?,?,?,?,?)',
                        (term['term_id'], term['evidence_id'], term.get('authority_native_code'), term.get('original_ko'), term.get('original_en'), canonical(term.get('aliases', [])), term['entity_scope'], term['scope_review_state'], canonical(term)))
        for binding in seed.get('compound_bindings', []):
            if not con.execute('SELECT 1 FROM compound WHERE full_inchikey=?', (binding['full_inchikey'],)).fetchone(): continue
            con.execute('INSERT OR IGNORE INTO compound_term_binding VALUES(?,?,?,?,?,?,?)',
                        tuple(binding.get(k) for k in ['full_inchikey', 'term_id', 'match_grade', 'review_state', 'matched_full_inchikey', 'identity_evidence_id', 'rationale']))

def report(con):
    def rows(sql):
        cursor = con.execute(sql); cols = [d[0] for d in cursor.description]
        return [dict(zip(cols, row)) for row in cursor.fetchall()]
    counts = {name: con.execute('SELECT count(*) FROM ' + name).fetchone()[0] for name in
              ['source_snapshot','raw_record','compound','structure_representation','compound_alias','food','food_context','observation','official_term','compound_term_binding','food_term_binding']}
    return {'generated_utc': datetime.now(timezone.utc).isoformat(), 'counts': counts,
            'compound_language_coverage': rows('SELECT * FROM v_compound_language_coverage')[0],
            'strict_projection_language_coverage': rows("SELECT count(DISTINCT o.full_inchikey) denominator_compounds,count(DISTINCT CASE WHEN EXISTS(SELECT 1 FROM compound_alias a WHERE a.full_inchikey=o.full_inchikey AND a.language='en') THEN o.full_inchikey END) source_en_present,count(DISTINCT CASE WHEN EXISTS(SELECT 1 FROM compound_alias a WHERE a.full_inchikey=o.full_inchikey AND a.language='ko') THEN o.full_inchikey END) source_ko_present,count(DISTINCT CASE WHEN b.review_state='verified' AND b.match_grade='exact_form' THEN o.full_inchikey END) authority_verified_exact_bilingual,count(DISTINCT CASE WHEN b.review_state IN('candidate','verified') AND length(trim(coalesce(t.original_ko,'')))>0 THEN o.full_inchikey END) authority_ko_candidate_or_verified FROM observation o LEFT JOIN compound_term_binding b USING(full_inchikey) LEFT JOIN official_term t USING(term_id) WHERE quality_state='accepted_strict'")[0],
            'observation_status_counts': rows('SELECT quality_state,count(*) row_count FROM observation GROUP BY quality_state'),
            'held_status_counts': rows('SELECT * FROM v_held_status_counts ORDER BY row_count DESC'),
            'identifier_conflicts': rows('SELECT * FROM v_identifier_conflicts'),
            'food_context_name_counts': rows('SELECT count(*) context_denominator,count(*) FILTER (WHERE original_common_name IS NOT NULL AND trim(original_common_name)!=\'\') source_name_populated,count(*) FILTER (WHERE original_part IS NOT NULL AND trim(original_part)!=\'\') part_populated,count(*) FILTER (WHERE taxon_state=\'authority_verified\') verified_taxon_contexts,count(*) FILTER (WHERE preparation_state=\'authority_verified\') verified_preparation_contexts FROM food_context')[0],
            'foreign_key_check': rows('PRAGMA foreign_key_check'), 'integrity_check': con.execute('PRAGMA integrity_check').fetchone()[0],
            'scope': 'Every staged audit row is retained in an isolated source projection. The separate 67,706-entry integrated reference population is not available in these inputs; its coverage is unknown. This database has not been promoted to an operational deployment.',
            'coverage_definition': 'Source language presence uses NFC text heuristics and is not canonical or authority verification. Exact bilingual coverage requires both nonblank primary-source labels, reviewed scope and exact form/context identity evidence.'}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--database', type=Path, default=HERE / 'food_quality.sqlite')
    parser.add_argument('--projection', type=Path, default=HERE / 'immutable_food_projection.ndjson.gz')
    parser.add_argument('--expected-sha256', required=True)
    parser.add_argument('--term-seed', type=Path)
    parser.add_argument('--report', type=Path, default=HERE / 'quality_report.json')
    args = parser.parse_args()
    if args.database.resolve().parent != HERE: raise ValueError('database must be in this owned isolated directory')
    con = connect(args.database)
    try:
        initialize(con)
        imported = import_projection(con, args.projection, args.expected_sha256)
        if args.term_seed: add_term_seed(con, json.loads(args.term_seed.read_bytes()))
        result = report(con); result['migration'] = imported
        args.report.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
        print(json.dumps({'migration': imported, 'counts': result['counts'], 'language': result['compound_language_coverage'], 'integrity': result['integrity_check']}, ensure_ascii=False))
    finally: con.close()

if __name__ == '__main__': main()
