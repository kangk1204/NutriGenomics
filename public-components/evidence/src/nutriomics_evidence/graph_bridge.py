"""Connect frozen food chemistry and regulatory documents without causal promotion."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from .graph import canonical, connect, edge, node, source, validate
from .sources import sha256


def _record(path, label, url):
    digest = sha256(path)
    return {'id': label + ':' + digest[:16], 'sha256': digest,
            'file': path.name, 'source_url': url, 'version': '2026-10-01 research snapshot',
            'redistribution_verified': False}


def integrate(database: Path, bridge: Path, regulatory: Path):
    db = connect(database)
    try:
        with db:
            composition = bridge / 'food_composition.csv.gz'
            identities = bridge / 'verified_compounds.csv'
            links = bridge / 'ctd_identity_links.csv'
            manifest = json.loads((bridge / 'bridge_source_manifest.json').read_text())
            food_source = _record(composition, 'USDA-Flavonoids-R03-3', manifest['sources'][0]['url'])
            food_source['upstream_provenance'] = manifest['sources']
            identity_source = _record(identities, 'PubChem-food-panel', 'https://pubchem.ncbi.nlm.nih.gov/')
            source(db, food_source); source(db, identity_source)
            compounds = {}
            with identities.open(newline='', encoding='utf-8') as stream:
                for row in csv.DictReader(stream):
                    identifier = 'InChIKey:' + row['compound_id']
                    compounds[row['compound_id']] = identifier
                    node(db, identifier, 'chemical', row['name'], {'full_inchikey': row['compound_id'], 'smiles': row['smiles'], 'structure_verified_by': 'PubChem/RDKit food-panel audit'})
                    for namespace, value in [('InChIKey', row['compound_id']), ('PubChem', row['pubchem_cid'])]:
                        db.execute('INSERT OR IGNORE INTO identifier_mapping VALUES (?,?,?,?,?,?)', (identifier, namespace, value, 'source_reported', identity_source['id'], 'Official PubChem record; full structure retained'))
            count = 0
            with gzip.open(composition, 'rt', newline='', encoding='utf-8') as stream:
                for row in csv.DictReader(stream):
                    compound = compounds[row['compound_id']]
                    food = 'USDA-Flav-R03-3:' + row['NDB_No']
                    node(db, food, 'food', row['Long_Desc'], {'active_source': food_source['id'], 'preparation': row['Long_Desc'], 'scientific_name': row['SciName'], 'cross_database_food_identity_verified': False})
                    amount = float(row['Flav_Val'])
                    if not math.isfinite(amount) or amount < 0:
                        raise ValueError('Invalid USDA flavonoid quantity')
                    # Hydrolyzed/aglycone-equivalent values remain analytical equivalents.
                    context = {'source': 'USDA Flavonoids R03-3', 'food_basis': row['food_basis'], 'amount': amount,
                               'unit': 'mg', 'basis': 'source_100g', 'source_quality_code': row['CC'],
                               'oral_dose_to_target': False, 'plasma_concentration': None}
                    identifier = hashlib.sha256(canonical([food, compound, row]).encode()).hexdigest()
                    status = 'measured_zero' if amount == 0 else 'measured'
                    db.execute('INSERT OR IGNORE INTO measurement VALUES (?,?,?,?,?,?,?,?,?)', (identifier, food, compound, amount, 'mg', 'source_100g', status, food_source['id'], canonical(row)))
                    edge(db, food, compound, 'contains_analytical_flavonoid', 'observed_composition', None, context, food_source['id'])
                    count += 1
            link_source = _record(links, 'CTD-PubChem-food-crossrefs', 'https://ctdbase.org/downloads/')
            source(db, link_source)
            link_count = 0
            with links.open(newline='', encoding='utf-8') as stream:
                for row in csv.DictReader(stream):
                    ctd_id = row['ctd_chemical_id']
                    node(db, ctd_id, 'chemical', row['ctd_name'])
                    # A CTD CID without an InChIKey cannot certify molecular-form equivalence.
                    edge(db, compounds[row['compound_id']], ctd_id, 'source_identifier_cross_reference', 'curated', None,
                         {'match_basis': row['match_basis'], 'ctd_full_structure_available': False, 'exact_structure_identity_claimed': False}, link_source['id'])
                    link_count += 1
            document = regulatory / 'regulatory_annotations.json'
            manifest_reg = json.loads((regulatory / 'source_manifest.json').read_text())
            if sha256(document) != manifest_reg['annotations_sha256']:
                raise ValueError('Regulatory annotations hash mismatch')
            payload = json.loads(document.read_text())
            db.execute('CREATE TABLE IF NOT EXISTS regulatory_annotation(id TEXT PRIMARY KEY, tier TEXT NOT NULL CHECK(tier=\'regulatory_annotation\'), metadata_json TEXT NOT NULL, source_id TEXT NOT NULL REFERENCES source(id))')
            reg_source = _record(document, 'FDA-eCFR-conditions', 'https://www.fda.gov/food/')
            source(db, reg_source)
            for row in payload['annotations']:
                if not row['source_checked'] or not row['all_configured_sources_verified']:
                    raise ValueError('Unverified regulatory source')
                for original in row['frozen_sources']:
                    file = regulatory / original['file']
                    if file.parent.resolve() != regulatory.resolve() or sha256(file) != original['sha256']:
                        raise ValueError('Regulatory original file/hash invalid')
                if row['tier'] != 'regulatory_annotation' or row['binding_affinity'] is not None or row['observed_individual_health_effect'] is not None:
                    raise ValueError('Regulatory conditions promoted to observed efficacy')
                db.execute('INSERT OR REPLACE INTO regulatory_annotation VALUES (?,?,?,?)', (row['id'], row['tier'], canonical(row), reg_source['id']))
            result = validate(db)
            if result['errors']:
                raise ValueError(result['errors'])
            result.update({'food_panel': {'compounds': len(compounds), 'composition_rows': count, 'ctd_identifier_crossrefs': link_count},
                           'regulatory_annotations': db.execute('SELECT count(*) FROM regulatory_annotation').fetchone()[0],
                           'cross_database_food_deduplication': 'Source-specific food IDs retained; no automatic name merge',
                           'FDA_conditions_are_efficacy_measurements': False})
        return result
    finally:
        db.close()
