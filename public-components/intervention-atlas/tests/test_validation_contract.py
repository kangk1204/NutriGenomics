"""Exact-function AST regression tests; no full scientific package import.

The exact validator and unchanged BH functions come from this component.
Native provenance I/O is loaded; SciPy and analysis routines are not imported.
"""
import ast
import hashlib
import json
import importlib.util
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'src/nutriomics_atlas/validation.py'
IO_SOURCE = ROOT / 'src/nutriomics_atlas/provenance.py'
spec = importlib.util.spec_from_file_location('validation_test_native_provenance', IO_SOURCE)
provenance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)
sha256 = provenance.sha256


namespace = dict(np=np, pd=pd, json=json, re=re, Path=Path, sha256=sha256)
for path, function in [(ROOT / 'src/nutriomics_atlas/metabolomics.py', 'bh'), (SOURCE, 'validate_output')]:
    node = next(n for n in ast.parse(path.read_bytes().decode('utf-8')).body
                if isinstance(n, ast.FunctionDef) and n.name == function)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
validate_output = namespace['validate_output']


@pytest.fixture
def table():
    return pd.DataFrame(dict(feature_id=['synthetic:one', 'synthetic:two'],
                             contrast=['arm_post_vs_pre'] * 2,
                             effect_estimate=[0.5, -0.2], standard_error=[0.1, 0.1],
                             confidence_low=[0.3, -0.4], confidence_high=[0.7, 0.0],
                             p_value=[0.01, 0.2], q_value=[0.02, 0.2],
                             q_value_study_family=[0.02, 0.2],
                             status=['estimable', 'estimable'], n_people=[3, 3]))


def write_fixture(directory, table, metadata=None):
    directory.mkdir(parents=True, exist_ok=True)
    table.to_csv(directory / 'metabolite_effects.tsv.gz', sep='\t', index=False, compression='gzip')
    record = dict(accession='ST001257', independent_people=19)
    if metadata:
        record.update(metadata)
    (directory / 'analysis_metadata.json').write_text(json.dumps(record), encoding='utf-8')
    return directory


def record(path):
    return {'path': str(path.resolve()), 'sha256': sha256(path), 'bytes': path.stat().st_size}


def test_valid_structural_table(tmp_path, table):
    assert validate_output(write_fixture(tmp_path, table))['errors'] == []


@pytest.mark.parametrize(('column', 'values'), [
    ('q_value_study_family', [0.0, 0.0]),
    ('q_value_study_family', [-3.0, 8.0]),
    ('q_value_study_family', [np.nan, 0.2]),
    ('q_value', [0.01, 0.2]),
    ('effect_estimate', [np.nan, -0.2]),
    ('effect_estimate', [np.inf, -0.2]),
    ('standard_error', [-0.1, 0.1]),
    ('standard_error', [np.nan, 0.1]),
    ('standard_error', [0.0, 0.1]),
    ('confidence_low', [0.9, -0.4]),
    ('confidence_high', [0.4, 0.0]),
    ('confidence_high', [np.nan, 0.0]),
    ('n_people', [0, 3]),
    ('n_people', [-1, 3]),
    ('n_people', [np.nan, 3]),
    ('n_people', [0.5, 3]),
    ('feature_id', ['', 'synthetic:two']),
    ('contrast', [' ', 'arm_post_vs_pre']),
    ('p_value', [-0.01, 0.2]),
    ('p_value', ['bad', 0.2]),
    ('p_value', [np.inf, 0.2]),
])
def test_malformed_estimable_rows_return_errors(tmp_path, table, column, values):
    table[column] = values
    result = validate_output(write_fixture(tmp_path, table))
    assert result['errors']
    assert result['patient_multiomics_fusion_confirmed'] is False


def test_header_only_rejected_without_exception(tmp_path, table):
    assert validate_output(write_fixture(tmp_path, table.iloc[:0]))['errors']


@pytest.mark.parametrize('zero_pairs', [False, True])
def test_explicit_nonestimable_rows_remain_supported(tmp_path, table, zero_pairs):
    columns = ['p_value', 'q_value', 'q_value_study_family', 'confidence_low', 'confidence_high']
    if zero_pairs:
        columns += ['effect_estimate', 'standard_error']
    for column in columns:
        table.loc[0, column] = np.nan
    table.loc[0, 'standard_error'] = np.nan if zero_pairs else 0
    table.loc[0, 'n_people'] = 0 if zero_pairs else 3
    table.loc[0, 'status'] = 'insufficient_pairs_or_constant_delta'
    result = validate_output(write_fixture(tmp_path, table))
    assert result['errors'] == []
    assert result['tables']['metabolite_effects.tsv.gz']['estimable_tests'] == 1


def test_missing_p_needs_explicit_nonestimable_status(tmp_path, table):
    for column in ['p_value', 'q_value', 'q_value_study_family']:
        table.loc[0, column] = np.nan
    assert validate_output(write_fixture(tmp_path, table))['errors']


def test_separate_contrast_q_retained_without_changing_study_family(tmp_path, table):
    table['contrast'] = ['arm_one', 'arm_two']
    table['q_value_contrast'] = [0.01, 0.2]
    assert validate_output(write_fixture(tmp_path, table))['errors'] == []


def test_negative_se_rejected_even_on_nonestimable_row(tmp_path, table):
    for column in ['p_value', 'q_value', 'q_value_study_family', 'confidence_low', 'confidence_high']:
        table.loc[0, column] = np.nan
    table.loc[0, 'status'] = 'insufficient_pairs_or_constant_delta'
    table.loc[0, 'standard_error'] = -0.1
    assert validate_output(write_fixture(tmp_path, table))['errors']


@pytest.mark.parametrize('samples', [0, -1, 2, 6.5, np.nan])
def test_invalid_gene_sample_counts(tmp_path, table, samples):
    table['n_samples'] = samples
    table['status'] = 'measured'
    for name in ['gene_effects.tsv.gz', 'pathway_effects.tsv.gz']:
        table.to_csv(tmp_path / name, sep='\t', index=False, compression='gzip')
    (tmp_path / 'analysis_metadata.json').write_text(json.dumps({'accession': 'GSE27385', 'independent_people': 3}))
    assert validate_output(tmp_path)['errors']


def test_native_filename_result_inventory_checked(tmp_path, table):
    directory = write_fixture(tmp_path, table, {'result_files': ['metabolite_effects.tsv.gz']})
    assert validate_output(directory)['errors'] == []
    write_fixture(tmp_path, table, {'result_files': ['metabolite_effects.tsv.gz', 'missing.tsv.gz']})
    assert validate_output(directory)['errors']


def test_native_hashed_result_inventory_checked(tmp_path, table):
    directory = write_fixture(tmp_path, table)
    receipt = record(directory / 'metabolite_effects.tsv.gz')
    metadata = json.loads((directory / 'analysis_metadata.json').read_text())
    metadata['result_files'] = [receipt]
    (directory / 'analysis_metadata.json').write_text(json.dumps(metadata))
    assert validate_output(directory)['errors'] == []
    receipt['sha256'] = '0' * 64
    (directory / 'analysis_metadata.json').write_text(json.dumps(metadata))
    assert validate_output(directory)['errors']


@pytest.mark.parametrize('records', [[], ['missing.tsv.gz'], [{'path': 'metabolite_effects.tsv.gz'}]])
def test_declared_result_inventory_cannot_omit_required_tables(tmp_path, table, records):
    assert validate_output(write_fixture(tmp_path, table, {'result_files': records}))['errors']


def test_hashed_result_inventory_covers_primary_table(tmp_path, table):
    directory = write_fixture(tmp_path, table)
    (directory / 'auxiliary.tsv').write_text('synthetic\n')
    metadata = json.loads((directory / 'analysis_metadata.json').read_text())
    metadata['result_files'] = ['metabolite_effects.tsv.gz', record(directory / 'auxiliary.tsv')]
    (directory / 'analysis_metadata.json').write_text(json.dumps(metadata))
    assert validate_output(directory)['errors']


def test_complete_native_run_receipt(tmp_path, table):
    directory = write_fixture(tmp_path, table)
    files = [directory / 'metabolite_effects.tsv.gz', directory / 'analysis_metadata.json']
    receipt = {'status': 'completed', 'returncode': 0, 'outputs': [record(path) for path in files]}
    (directory / 'run_manifest.json').write_text(json.dumps(receipt))
    assert validate_output(directory)['errors'] == []
    receipt['outputs'][0]['path'] = r'C:\synthetic\metabolite_effects.tsv.gz'
    (directory / 'run_manifest.json').write_text(json.dumps(receipt))
    assert validate_output(directory)['errors'] == []


@pytest.mark.parametrize('missing_metadata', [False, True])
def test_run_receipt_cannot_omit_required_hashes(tmp_path, table, missing_metadata):
    directory = write_fixture(tmp_path, table)
    files = [directory / 'metabolite_effects.tsv.gz'] if missing_metadata else []
    receipt = {'status': 'completed', 'returncode': 0, 'outputs': [record(path) for path in files]}
    (directory / 'run_manifest.json').write_text(json.dumps(receipt))
    assert validate_output(directory)['errors']


def test_run_receipt_checks_every_declared_record(tmp_path, table):
    directory = write_fixture(tmp_path, table)
    receipt = {'status': 'completed', 'returncode': 0,
               'outputs': [record(directory / 'metabolite_effects.tsv.gz'), record(directory / 'analysis_metadata.json'),
                           {'path': 'missing_aux.tsv', 'sha256': '0' * 64}]}
    (directory / 'run_manifest.json').write_text(json.dumps(receipt))
    assert validate_output(directory)['errors']


def test_native_pathway_nonestimable_row_needs_no_gene_ci(tmp_path, table):
    table['n_samples'] = 6
    table['status'] = 'measured'
    table.to_csv(tmp_path / 'gene_effects.tsv.gz', sep='\t', index=False, compression='gzip')
    pathway = table.drop(columns=['confidence_low', 'confidence_high']).copy()
    for column in ['p_value', 'q_value', 'q_value_study_family', 'effect_estimate', 'standard_error']:
        pathway.loc[0, column] = np.nan
    pathway.loc[0, 'status'] = 'not_estimable'
    pathway.to_csv(tmp_path / 'pathway_effects.tsv.gz', sep='\t', index=False, compression='gzip')
    (tmp_path / 'analysis_metadata.json').write_text(json.dumps({'accession': 'GSE27385', 'independent_people': 3}))
    result = validate_output(tmp_path)
    assert result['errors'] == []
    assert result['tables']['pathway_effects.tsv.gz']['estimable_tests'] == 1


@pytest.mark.parametrize('cohort', [0, -1, 2.5, None, 'bad', True])
def test_invalid_cohort_count_returns_error(tmp_path, table, cohort):
    assert validate_output(write_fixture(tmp_path, table, {'independent_people': cohort}))['errors']
# Append these two controls to the existing validation contract test module.

def test_r_measured_row_missing_p_is_reported(tmp_path, table):
    # Native R gene exporters label every row measured. A missing test must
    # be reported as inconsistent rather than silently relabeled or filled.
    table['status'] = 'measured'
    for column in ['p_value', 'q_value', 'q_value_study_family']:
        table.loc[0, column] = np.nan
    assert validate_output(write_fixture(tmp_path, table))['errors']


def test_all_explicit_nonestimable_rows_preserved(tmp_path, table):
    # An all-missing feature has zero complete pairs, not zero-valued effects.
    for column in ['effect_estimate', 'standard_error', 'p_value', 'q_value',
                   'q_value_study_family', 'confidence_low', 'confidence_high']:
        table[column] = np.nan
    table['n_people'] = 0
    table['status'] = 'insufficient_pairs_or_constant_delta'
    result = validate_output(write_fixture(tmp_path, table))
    assert result['errors'] == []
    summary = result['tables']['metabolite_effects.tsv.gz']
    assert summary['rows'] == 2 and summary['estimable_tests'] == 0
    assert summary['q_below_005'] == 0
