import json
from pathlib import Path
import pytest
from nutriomics_evidence.graph import connect
from nutriomics_evidence.graph_bridge import integrate

ROOT = Path(__file__).resolve().parents[2]

def test_frozen_food_regulatory_tiers(tmp_path):
    bridge = ROOT / 'nutriomics-compound-target/data/food_bridge'
    regulatory = Path(__file__).resolve().parents[1] / 'data/fda_claims'
    if not bridge.exists() or not regulatory.exists():
        pytest.skip('Frozen private research inputs are not distributed in Git')
    database = tmp_path / 'graph.sqlite'
    result = integrate(database, bridge, regulatory)
    assert result['food_panel'] == {'compounds': 12, 'composition_rows': 2628, 'ctd_identifier_crossrefs': 4}
    assert result['regulatory_annotations'] == 4
    db = connect(database, readonly=True)
    assert not db.execute("SELECT * FROM edge WHERE relation IN ('direct_binding','dti_positive')").fetchall()
    assert not db.execute("SELECT * FROM identifier_mapping WHERE relation='exact_structure'").fetchall()
    assert all(json.loads(row[0])['observed_individual_health_effect'] is None for row in db.execute('SELECT metadata_json FROM regulatory_annotation'))
    db.close()
    # The same frozen input may be rerun without duplicating measurements/edges.
    repeat = integrate(database, bridge, regulatory)
    assert repeat['tables'] == result['tables']
