import bz2
import csv
import hashlib
import json
from pathlib import Path
import pytest
from nutriomics_evidence.food_compound_extension import verify_structure,curated_parts,stage,strict_analyte_status,normalized_analyte_name

def test_mislabeled_payload_binding_exact_key_and_mixture_guard():
    row=('7','FDB000007','Ethanol','CCO','InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3','LFQSCWFLJHTTHZ-UHFFFAOYSA-N')
    assert verify_structure(row)['standard_structure_verified']
    bad=verify_structure(row[:-1]+('AAAAAAAAAAAAAA-UHFFFAOYSA-N',))
    assert not bad['standard_structure_verified']
    mixed=verify_structure(row[:3]+('CCO.O',)+row[4:])
    assert not mixed['single_component'] and not mixed['standard_structure_verified']
    aggregate=verify_structure(row[:2]+('Total carbohydrates',)+row[3:])
    assert aggregate['aggregate_concept'] and not aggregate['standard_structure_verified']

def test_published_edible_parts_remain_pairs_not_global_part_allowlist(tmp_path):
    path=tmp_path/'parts.Rmd';path.write_text('desired_orig_food_parts = data.frame(food_id = c(6, 9), desired_part = c("Bulb", "Leaf"))')
    rows=curated_parts(path)
    assert list(rows.itertuples(index=False,name=None))==[('6','Bulb'),('9','Leaf')]
    assert ('6','Leaf') not in set(rows.itertuples(index=False,name=None))
    path.write_text('desired_orig_food_parts = data.frame(food_id = c(6, 9), desired_part = c("Bulb"))')
    with pytest.raises(ValueError,match='unpaired'):curated_parts(path)


def test_actual_stage_excludes_predictions_wrong_parts_and_impossible_values(tmp_path):
    raw=tmp_path/'raw';raw.mkdir()
    def write(name,rows):
        target=raw/name
        opener=bz2.open if name.endswith('.bz2') else open
        with opener(target,'wt',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        return {'name':name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    compounds=[dict(id='7',public_id='FDB000007',name='Ethanol',cas_number='CCO',moldb_inchikey='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3',moldb_smiles='LFQSCWFLJHTTHZ-UHFFFAOYSA-N')]
    foods=[dict(id='6',name='Apple',food_group='Fruits',food_type='Type 1',export_to_foodb='1',public_id='FOOD00006')]
    base=dict(id='1',source_id='7',source_type='Compound',food_id='6',orig_food_common_name='Apple',orig_food_part='Fruit',orig_content='1',orig_min='',orig_max='',standard_content='',orig_unit='mg/100g',orig_citation='',citation='USDA',citation_type='DATABASE',orig_method='',preparation_type='Raw')
    rows=[base,
      dict(base,id='2',citation='HMDB',citation_type='PREDICTED'),
      dict(base,id='3',citation='DUKE',orig_food_part='Leaf'),
      dict(base,id='4',citation='DUKE',orig_content=''),
      dict(base,id='5',orig_content='0'),
      dict(base,id='6',source_type='Nutrient'),
      dict(base,id='7',citation='MANUAL',citation_type='ARTICLE'),
      dict(base,id='8',orig_content='100001'),
      dict(base,id='9',citation='PHYTOHUB',orig_content='')]
    sources=[write('Compound.csv.bz2',compounds),write('Content.csv.bz2',rows),write('Food.csv',foods)]
    parts=raw/'edible_part_rules.Rmd';parts.write_text('desired_orig_food_parts = data.frame(food_id = c(6), desired_part = c("Fruit"))')
    sources.append({'name':parts.name,'sha256':hashlib.sha256(parts.read_bytes()).hexdigest()})
    (tmp_path/'acquisition.json').write_text(json.dumps({'sources':sources}))
    receipt=stage(tmp_path)
    assert receipt['verified_food_single_compound_standard_count']==1
    assert receipt['eligible_counts']['food_compound_extension_eligible_occurrences']==2
    assert receipt['all_checks_zero']
    assert receipt['rejection_counts']['predicted_not_observed_food']==1
    assert receipt['rejection_counts']['duke_part_not_in_published_food_specific_curation']==1
    assert receipt['rejection_counts']['impossible_mass_per100g']==1
    assert not receipt['meets_2000_occurrence_verified_standard_target']


def test_strict_identity_distinguishes_quantitative_presence_and_assignments():
    row=dict(positive_value=1,standard_structure_verified=True,single_component=True,standard_unspecified_stereo=False,
        source_pmid=None,name='Ethanol',orig_source_id='ETOH',orig_source_name='ETHANOL',orig_unit='mg/100g',orig_method='',
        citation='DUKE',orig_citation=None,published_edible_part_verified=True,eligibility='eligible_positive_curated_edible_part')
    assert strict_analyte_status(row)=='strict_positive_source_named_single_analyte'
    assert strict_analyte_status(dict(row,source_pmid='30994344'))=='milk_assay_to_structure_mapping_unresolved'
    assert strict_analyte_status(dict(row,name='Total flavonoids'))=='generic_class_or_lipid_species_assignment'
    assert strict_analyte_status(dict(row,name='Vitamin D',orig_source_name='Vitamin D'))=='generic_class_or_lipid_species_assignment'
    assert strict_analyte_status(dict(row,name='Vitamin D3',orig_source_name='Vitamin D3'))=='strict_positive_source_named_single_analyte'
    assert strict_analyte_status(dict(row,standard_unspecified_stereo=True))=='standard_stereochemistry_unspecified'
    assert strict_analyte_status(dict(row,orig_source_name=None))=='original_assayed_analyte_id_or_name_missing'
    assert strict_analyte_status(dict(row,orig_source_name='Methanol'))=='assayed_analyte_name_not_exact_standard_identity'
    qualitative=dict(row,positive_value=None,orig_unit=None,eligibility='eligible_qualitative_curated_edible_part')
    assert strict_analyte_status(qualitative)=='qualitative_database_assignment_without_specific_source'
    assert strict_analyte_status(dict(qualitative,citation='12345678',citation_type='ARTICLE'))=='strict_source_reported_qualitative_single_analyte'
    assert normalized_analyte_name('(+)-Catechin')!=normalized_analyte_name('(-)-Catechin')
    assert strict_analyte_status(row,single_atom=True)=='element_total_not_resolved_molecular_form'
