PRAGMA foreign_keys=ON;
CREATE TABLE IF NOT EXISTS source_snapshot(
 snapshot_id TEXT PRIMARY KEY, protocol TEXT NOT NULL, frozen_utc TEXT NOT NULL,
 source_repository TEXT NOT NULL, source_commit TEXT NOT NULL, input_identity_json TEXT NOT NULL,
 projection_sha256 TEXT NOT NULL CHECK(length(projection_sha256)=64), import_complete INTEGER NOT NULL DEFAULT 0 CHECK(import_complete IN(0,1))
);
CREATE TABLE IF NOT EXISTS source_file(
 snapshot_id TEXT NOT NULL REFERENCES source_snapshot, source_table TEXT NOT NULL,
 sha256 TEXT NOT NULL CHECK(length(sha256)=64), bytes INTEGER NOT NULL CHECK(bytes>=0), row_count INTEGER NOT NULL CHECK(row_count>=0),
 PRIMARY KEY(snapshot_id,source_table)
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS raw_record(
 snapshot_id TEXT NOT NULL, source_table TEXT NOT NULL, source_row_id TEXT NOT NULL,
 raw_sha256 TEXT NOT NULL CHECK(length(raw_sha256)=64), raw_json TEXT NOT NULL,
 PRIMARY KEY(snapshot_id,source_table,source_row_id),
 FOREIGN KEY(snapshot_id,source_table) REFERENCES source_file
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS compound(
 full_inchikey TEXT PRIMARY KEY CHECK(length(full_inchikey)=27 AND substr(full_inchikey,15,1)='-' AND substr(full_inchikey,26,1)='-' AND full_inchikey NOT GLOB '*[^A-Z-]*'),
 source_structure_verified INTEGER NOT NULL CHECK(source_structure_verified IN(0,1)),
 source_single_component INTEGER NOT NULL CHECK(source_single_component IN(0,1)),
 source_stereo_unspecified INTEGER NOT NULL CHECK(source_stereo_unspecified IN(0,1))
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS structure_representation(
 representation_id TEXT PRIMARY KEY, full_inchikey TEXT NOT NULL REFERENCES compound,
 snapshot_id TEXT NOT NULL REFERENCES source_snapshot, original_inchi TEXT, original_isomeric_smiles TEXT
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS compound_identifier_assertion(
 snapshot_id TEXT NOT NULL REFERENCES source_snapshot, namespace TEXT NOT NULL, native_id TEXT NOT NULL,
 full_inchikey TEXT NOT NULL REFERENCES compound,
 PRIMARY KEY(snapshot_id,namespace,native_id,full_inchikey)
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS compound_alias(
 full_inchikey TEXT NOT NULL REFERENCES compound, original_name TEXT NOT NULL,
 search_name TEXT NOT NULL, language TEXT NOT NULL CHECK(language IN('en','ko','und')),
 snapshot_id TEXT NOT NULL REFERENCES source_snapshot, name_role TEXT NOT NULL DEFAULT 'source_alias' CHECK(name_role='source_alias'),
 PRIMARY KEY(full_inchikey,original_name,language,snapshot_id)
) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS compound_alias_lookup ON compound_alias(search_name);
CREATE TABLE IF NOT EXISTS food(
 snapshot_id TEXT NOT NULL REFERENCES source_snapshot, native_food_id TEXT NOT NULL,
 PRIMARY KEY(snapshot_id,native_food_id)
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS food_identifier_assertion(
 snapshot_id TEXT NOT NULL, native_food_id TEXT NOT NULL, namespace TEXT NOT NULL, public_id TEXT NOT NULL,
 PRIMARY KEY(snapshot_id,native_food_id,namespace,public_id),
 FOREIGN KEY(snapshot_id,native_food_id) REFERENCES food
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS food_context(
 context_id TEXT PRIMARY KEY, snapshot_id TEXT NOT NULL, native_food_id TEXT NOT NULL,
 original_common_name TEXT, original_part TEXT,
 taxon_id TEXT, preparation_code TEXT,
 taxon_state TEXT NOT NULL DEFAULT 'not_provided' CHECK(taxon_state IN('not_provided','source_asserted','authority_verified')),
 preparation_state TEXT NOT NULL DEFAULT 'not_provided' CHECK(preparation_state IN('not_provided','source_asserted','authority_verified')),
 FOREIGN KEY(snapshot_id,native_food_id) REFERENCES food
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS observation(
 snapshot_id TEXT NOT NULL, source_row_id TEXT NOT NULL, raw_table TEXT NOT NULL DEFAULT 'audit' CHECK(raw_table='audit'),
 full_inchikey TEXT NOT NULL REFERENCES compound, context_id TEXT NOT NULL REFERENCES food_context,
 positive_value REAL, original_unit TEXT, normalized_unit_code TEXT,
 original_method TEXT, original_analyte_id TEXT, original_analyte_name TEXT,
 source_citation TEXT, source_citation_type TEXT, source_pmid TEXT,
 strict_status TEXT NOT NULL, quality_state TEXT NOT NULL CHECK(quality_state IN('accepted_strict','retained_held')),
 PRIMARY KEY(snapshot_id,source_row_id),
 FOREIGN KEY(snapshot_id,raw_table,source_row_id) REFERENCES raw_record,
 CHECK(quality_state!='accepted_strict' OR (strict_status='strict_positive_source_named_single_analyte' AND positive_value IS NOT NULL AND typeof(positive_value) IN('integer','real') AND positive_value>0 AND positive_value<=1.7976931348623157e308 AND original_unit IS NOT NULL AND original_unit='mg/100g')),
 CHECK(normalized_unit_code IS NULL OR (original_unit='mg/100g' AND normalized_unit_code='mg_per_100g'))
) WITHOUT ROWID;
CREATE INDEX IF NOT EXISTS observation_compound ON observation(full_inchikey);
CREATE INDEX IF NOT EXISTS observation_context ON observation(context_id);
CREATE TRIGGER IF NOT EXISTS guard_strict_structure BEFORE INSERT ON observation
WHEN NEW.quality_state='accepted_strict' AND NOT EXISTS(
 SELECT 1 FROM compound WHERE full_inchikey=NEW.full_inchikey AND source_structure_verified=1 AND source_single_component=1 AND source_stereo_unspecified=0)
BEGIN SELECT RAISE(ABORT,'strict observation requires verified single-component specified structure'); END;
CREATE TABLE IF NOT EXISTS authority_evidence(
 evidence_id TEXT PRIMARY KEY, authority TEXT NOT NULL, source_uri TEXT NOT NULL,
 source_version TEXT NOT NULL, retrieved_utc TEXT NOT NULL,
 verification_method TEXT NOT NULL CHECK(verification_method IN('direct_primary_source','delegated_primary_source','unverified')),
 evidence_sha256 TEXT NOT NULL CHECK(length(evidence_sha256)=64), evidence_json TEXT NOT NULL
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS official_term(
 term_id TEXT PRIMARY KEY, evidence_id TEXT NOT NULL REFERENCES authority_evidence,
 authority_native_code TEXT, original_ko TEXT, original_en TEXT, aliases_json TEXT NOT NULL DEFAULT '[]',
 entity_scope TEXT NOT NULL CHECK(entity_scope IN('compound_concept','food_context','organism','mixture')),
 scope_review_state TEXT NOT NULL CHECK(scope_review_state IN('reviewed','pending','conflict')),
 raw_source_json TEXT NOT NULL
) WITHOUT ROWID;
CREATE TABLE IF NOT EXISTS compound_term_binding(
 full_inchikey TEXT NOT NULL REFERENCES compound, term_id TEXT NOT NULL REFERENCES official_term,
 match_grade TEXT NOT NULL CHECK(match_grade IN('exact_form','broader_concept','unresolved','conflict')),
 review_state TEXT NOT NULL CHECK(review_state IN('verified','candidate','blocked')),
 matched_full_inchikey TEXT, identity_evidence_id TEXT REFERENCES authority_evidence,
 rationale TEXT NOT NULL,
 PRIMARY KEY(full_inchikey,term_id),
 CHECK(review_state!='verified' OR (match_grade='exact_form' AND matched_full_inchikey=full_inchikey AND identity_evidence_id IS NOT NULL))
) WITHOUT ROWID;
CREATE TRIGGER IF NOT EXISTS guard_verified_compound_term BEFORE INSERT ON compound_term_binding
WHEN NEW.review_state='verified' AND NOT EXISTS(
 SELECT 1 FROM official_term t JOIN authority_evidence a ON a.evidence_id=t.evidence_id
 JOIN authority_evidence i ON i.evidence_id=NEW.identity_evidence_id
 WHERE t.term_id=NEW.term_id AND t.scope_review_state='reviewed' AND t.entity_scope='compound_concept'
 AND length(trim(coalesce(t.original_ko,'')))>0 AND length(trim(coalesce(t.original_en,'')))>0
 AND a.verification_method='direct_primary_source' AND i.verification_method='direct_primary_source'
 AND json_extract(i.evidence_json,'$.full_inchikey')=NEW.full_inchikey
 AND json_extract(i.evidence_json,'$.form_review_passed')=1)
BEGIN SELECT RAISE(ABORT,'verified bilingual binding requires reviewed primary-source terms and exact-form identity evidence'); END;
CREATE TABLE IF NOT EXISTS food_term_binding(
 context_id TEXT NOT NULL REFERENCES food_context, term_id TEXT NOT NULL REFERENCES official_term,
 match_grade TEXT NOT NULL CHECK(match_grade IN('exact_context','broader_organism','unresolved','conflict')),
 review_state TEXT NOT NULL CHECK(review_state IN('verified','candidate','blocked')),
 taxon_match INTEGER NOT NULL CHECK(taxon_match IN(0,1)), part_match INTEGER NOT NULL CHECK(part_match IN(0,1)), preparation_match INTEGER NOT NULL CHECK(preparation_match IN(0,1)),
 identity_evidence_id TEXT REFERENCES authority_evidence, rationale TEXT NOT NULL,
 PRIMARY KEY(context_id,term_id),
 CHECK(review_state!='verified' OR (match_grade='exact_context' AND taxon_match=1 AND part_match=1 AND preparation_match=1 AND identity_evidence_id IS NOT NULL))
) WITHOUT ROWID;
CREATE TRIGGER IF NOT EXISTS guard_verified_food_term BEFORE INSERT ON food_term_binding
WHEN NEW.review_state='verified' AND NOT EXISTS(
 SELECT 1 FROM official_term t JOIN authority_evidence a ON a.evidence_id=t.evidence_id
 JOIN authority_evidence i ON i.evidence_id=NEW.identity_evidence_id
 WHERE t.term_id=NEW.term_id AND t.entity_scope='food_context' AND t.scope_review_state='reviewed'
 AND a.verification_method='direct_primary_source'
 AND length(trim(coalesce(t.original_ko,'')))>0 AND length(trim(coalesce(t.original_en,'')))>0
 AND i.verification_method='direct_primary_source'
 AND json_extract(i.evidence_json,'$.food_context_id')=NEW.context_id
 AND json_extract(i.evidence_json,'$.taxon_match')=1
 AND json_extract(i.evidence_json,'$.part_match')=1
 AND json_extract(i.evidence_json,'$.preparation_match')=1)
BEGIN SELECT RAISE(ABORT,'official nonblank food labels alone are insufficient for exact context binding'); END;
CREATE VIEW IF NOT EXISTS v_compound_language_coverage AS
SELECT count(*) denominator_compounds,
 sum(EXISTS(SELECT 1 FROM compound_alias a WHERE a.full_inchikey=c.full_inchikey AND a.language='en')) source_en_present,
 sum(EXISTS(SELECT 1 FROM compound_alias a WHERE a.full_inchikey=c.full_inchikey AND a.language='ko')) source_ko_present,
 sum(EXISTS(SELECT 1 FROM compound_term_binding b JOIN official_term t USING(term_id) WHERE b.full_inchikey=c.full_inchikey AND b.review_state IN('candidate','verified') AND length(trim(coalesce(t.original_ko,'')))>0)) authority_ko_candidate_or_verified,
 sum(EXISTS(SELECT 1 FROM compound_term_binding b WHERE b.full_inchikey=c.full_inchikey AND b.review_state='verified' AND b.match_grade='exact_form')) authority_verified_exact_bilingual
FROM compound c;
CREATE VIEW IF NOT EXISTS v_held_status_counts AS SELECT strict_status,count(*) row_count FROM observation WHERE quality_state='retained_held' GROUP BY strict_status;
CREATE VIEW IF NOT EXISTS v_identifier_conflicts AS SELECT snapshot_id,namespace,native_id,count(DISTINCT full_inchikey) compound_keys FROM compound_identifier_assertion GROUP BY snapshot_id,namespace,native_id HAVING count(DISTINCT full_inchikey)>1;
CREATE TRIGGER IF NOT EXISTS raw_record_no_update BEFORE UPDATE ON raw_record
BEGIN SELECT RAISE(ABORT,'raw source records are immutable; register a new source snapshot'); END;
CREATE TRIGGER IF NOT EXISTS raw_record_no_delete BEFORE DELETE ON raw_record
BEGIN SELECT RAISE(ABORT,'raw source records are immutable; discard the isolated DB for rollback'); END;
CREATE TRIGGER IF NOT EXISTS guard_observation_source_insert BEFORE INSERT ON observation
WHEN NOT EXISTS(SELECT 1 FROM raw_record r JOIN food_context f ON f.context_id=NEW.context_id
 WHERE r.snapshot_id=NEW.snapshot_id AND r.source_table='audit' AND r.source_row_id=NEW.source_row_id
 AND json_extract(r.raw_json,'$.full_inchikey')=NEW.full_inchikey
 AND json_extract(r.raw_json,'$.food_id')=f.native_food_id AND f.snapshot_id=NEW.snapshot_id
 AND json_extract(r.raw_json,'$.orig_food_common_name') IS f.original_common_name
 AND json_extract(r.raw_json,'$.orig_food_part') IS f.original_part
 AND json_extract(r.raw_json,'$.positive_value') IS NEW.positive_value
 AND CAST(json_extract(r.raw_json,'$.orig_unit') AS TEXT) IS NEW.original_unit
 AND CAST(json_extract(r.raw_json,'$.orig_method') AS TEXT) IS NEW.original_method
 AND CAST(json_extract(r.raw_json,'$.orig_source_id') AS TEXT) IS NEW.original_analyte_id
 AND CAST(json_extract(r.raw_json,'$.orig_source_name') AS TEXT) IS NEW.original_analyte_name
 AND CAST(json_extract(r.raw_json,'$.citation') AS TEXT) IS NEW.source_citation
 AND CAST(json_extract(r.raw_json,'$.citation_type') AS TEXT) IS NEW.source_citation_type
 AND CAST(json_extract(r.raw_json,'$.source_pmid') AS TEXT) IS NEW.source_pmid
 AND json_extract(r.raw_json,'$.strict_status')=NEW.strict_status
 AND NEW.quality_state=CASE WHEN json_extract(r.raw_json,'$.strict_status')='strict_positive_source_named_single_analyte' THEN 'accepted_strict' ELSE 'retained_held' END
 AND (NEW.quality_state!='accepted_strict' OR (
   json_extract(r.raw_json,'$.standard_structure_verified')=1
   AND json_extract(r.raw_json,'$.single_component')=1
   AND json_extract(r.raw_json,'$.standard_unspecified_stereo')=0)))
BEGIN SELECT RAISE(ABORT,'observation must preserve source identity, context, measurements, provenance and eligibility'); END;
CREATE TRIGGER IF NOT EXISTS guard_observation_source_update BEFORE UPDATE ON observation
WHEN NOT EXISTS(SELECT 1 FROM raw_record r JOIN food_context f ON f.context_id=NEW.context_id
 WHERE r.snapshot_id=NEW.snapshot_id AND r.source_table='audit' AND r.source_row_id=NEW.source_row_id
 AND json_extract(r.raw_json,'$.full_inchikey')=NEW.full_inchikey
 AND json_extract(r.raw_json,'$.food_id')=f.native_food_id AND f.snapshot_id=NEW.snapshot_id
 AND json_extract(r.raw_json,'$.orig_food_common_name') IS f.original_common_name
 AND json_extract(r.raw_json,'$.orig_food_part') IS f.original_part
 AND json_extract(r.raw_json,'$.positive_value') IS NEW.positive_value
 AND CAST(json_extract(r.raw_json,'$.orig_unit') AS TEXT) IS NEW.original_unit
 AND CAST(json_extract(r.raw_json,'$.orig_method') AS TEXT) IS NEW.original_method
 AND CAST(json_extract(r.raw_json,'$.orig_source_id') AS TEXT) IS NEW.original_analyte_id
 AND CAST(json_extract(r.raw_json,'$.orig_source_name') AS TEXT) IS NEW.original_analyte_name
 AND CAST(json_extract(r.raw_json,'$.citation') AS TEXT) IS NEW.source_citation
 AND CAST(json_extract(r.raw_json,'$.citation_type') AS TEXT) IS NEW.source_citation_type
 AND CAST(json_extract(r.raw_json,'$.source_pmid') AS TEXT) IS NEW.source_pmid
 AND json_extract(r.raw_json,'$.strict_status')=NEW.strict_status
 AND NEW.quality_state=CASE WHEN json_extract(r.raw_json,'$.strict_status')='strict_positive_source_named_single_analyte' THEN 'accepted_strict' ELSE 'retained_held' END
 AND (NEW.quality_state!='accepted_strict' OR (
   json_extract(r.raw_json,'$.standard_structure_verified')=1
   AND json_extract(r.raw_json,'$.single_component')=1
   AND json_extract(r.raw_json,'$.standard_unspecified_stereo')=0)))
BEGIN SELECT RAISE(ABORT,'observation must preserve source identity, context, measurements, provenance and eligibility'); END;
CREATE TRIGGER IF NOT EXISTS guard_verified_compound_term_update BEFORE UPDATE ON compound_term_binding
WHEN NEW.review_state='verified' AND NOT EXISTS(
 SELECT 1 FROM official_term t JOIN authority_evidence a ON a.evidence_id=t.evidence_id
 JOIN authority_evidence i ON i.evidence_id=NEW.identity_evidence_id
 WHERE t.term_id=NEW.term_id AND t.scope_review_state='reviewed' AND t.entity_scope='compound_concept'
 AND length(trim(coalesce(t.original_ko,'')))>0 AND length(trim(coalesce(t.original_en,'')))>0
 AND a.verification_method='direct_primary_source' AND i.verification_method='direct_primary_source'
 AND json_extract(i.evidence_json,'$.full_inchikey')=NEW.full_inchikey
 AND json_extract(i.evidence_json,'$.form_review_passed')=1)
BEGIN SELECT RAISE(ABORT,'verified bilingual update requires exact primary-source form evidence'); END;
CREATE TRIGGER IF NOT EXISTS guard_verified_food_term_update BEFORE UPDATE ON food_term_binding
WHEN NEW.review_state='verified' AND NOT EXISTS(
 SELECT 1 FROM official_term t JOIN authority_evidence a ON a.evidence_id=t.evidence_id
 JOIN authority_evidence i ON i.evidence_id=NEW.identity_evidence_id
 WHERE t.term_id=NEW.term_id AND t.entity_scope='food_context' AND t.scope_review_state='reviewed'
 AND a.verification_method='direct_primary_source'
 AND length(trim(coalesce(t.original_ko,'')))>0 AND length(trim(coalesce(t.original_en,'')))>0
 AND i.verification_method='direct_primary_source'
 AND json_extract(i.evidence_json,'$.food_context_id')=NEW.context_id
 AND json_extract(i.evidence_json,'$.taxon_match')=1
 AND json_extract(i.evidence_json,'$.part_match')=1
 AND json_extract(i.evidence_json,'$.preparation_match')=1)
BEGIN SELECT RAISE(ABORT,'verified food update requires matching taxon part and preparation evidence'); END;

CREATE TRIGGER IF NOT EXISTS authority_evidence_no_update BEFORE UPDATE ON authority_evidence
BEGIN SELECT RAISE(ABORT,'authority evidence is immutable; register a new evidence version'); END;
CREATE TRIGGER IF NOT EXISTS authority_evidence_no_delete BEFORE DELETE ON authority_evidence
BEGIN SELECT RAISE(ABORT,'authority evidence is immutable; discard the isolated DB for rollback'); END;
CREATE TRIGGER IF NOT EXISTS official_term_no_update BEFORE UPDATE ON official_term
BEGIN SELECT RAISE(ABORT,'official source terms are immutable; register a new term version'); END;
CREATE TRIGGER IF NOT EXISTS official_term_no_delete BEFORE DELETE ON official_term
BEGIN SELECT RAISE(ABORT,'official source terms are immutable; discard the isolated DB for rollback'); END;
CREATE TRIGGER IF NOT EXISTS food_context_source_no_update BEFORE UPDATE ON food_context
WHEN NEW.context_id!=OLD.context_id OR NEW.snapshot_id!=OLD.snapshot_id OR NEW.native_food_id!=OLD.native_food_id
 OR NEW.original_common_name IS NOT OLD.original_common_name OR NEW.original_part IS NOT OLD.original_part
BEGIN SELECT RAISE(ABORT,'original food identity, name and part are immutable source assertions'); END;
