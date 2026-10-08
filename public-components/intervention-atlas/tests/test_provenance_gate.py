"""Tiny receipt/status/CLI probes of actual applied sources, isolated by AST.

These tests do not import the SciPy-dependent analysis package or run inference.
Native functions and I/O helpers are loaded from this component only.
"""
import argparse
import ast
import importlib.util
import json
import re
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'src/nutriomics_atlas'
spec = importlib.util.spec_from_file_location('native_atlas_provenance_probe', SRC / 'provenance.py')
provenance = importlib.util.module_from_spec(spec)
spec.loader.exec_module(provenance)


def exact_function(path, name, namespace):
    node = next(n for n in ast.parse(path.read_bytes().decode('utf-8')).body
                if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace[name]


namespace = dict(np=np, pd=pd, json=json, re=re, Path=Path, sha256=provenance.sha256)
exact_function(SRC / 'metabolomics.py', 'bh', namespace)
validate_output = exact_function(SRC / 'validation.py', 'validate_output', namespace)
cli_namespace = dict(argparse=argparse, json=json, Path=Path, __package__='nutriomics_atlas')
cli_main = exact_function(SRC / 'cli.py', 'main', cli_namespace)


def prepare(directory, accession='ST001257', version=True):
    directory.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(dict(feature_id=['synthetic:one', 'synthetic:two'],
                             contrast=['arm_post_vs_pre'] * 2,
                             effect_estimate=[0.5, -0.2], standard_error=[0.1, 0.1],
                             confidence_low=[0.3, -0.4], confidence_high=[0.7, 0.0],
                             p_value=[0.01, 0.2], q_value=[0.02, 0.2],
                             q_value_study_family=[0.02, 0.2],
                             status=['estimable', 'estimable'], n_people=[3, 3]))
    metadata = {'accession': accession, 'independent_people': 3}
    if accession == 'ST001257':
        filenames = ['metabolite_effects.tsv.gz']
        if version:
            metadata['analysis_version'] = 'st001257-paired-processed-v1'
    else:
        filenames = ['gene_effects.tsv.gz', 'pathway_effects.tsv.gz']
        table['n_samples'] = 6
        table['status'] = 'measured'
        if version:
            metadata['analysis_version'] = 'gse127530-limma-voom-final-v1'
    for name in filenames:
        table.to_csv(directory / name, sep='\t', index=False, compression='gzip')
    save_metadata(directory, metadata)
    return metadata, filenames


def save_metadata(directory, metadata):
    (directory / 'analysis_metadata.json').write_text(json.dumps(metadata), encoding='utf-8')


def add_metadata_hashes(directory, metadata, filenames):
    metadata['result_files'] = [provenance.file_record(directory / name) for name in filenames]
    save_metadata(directory, metadata)


def add_run_hashes(directory, filenames, status='completed', returncode=0):
    receipt = {'status': status, 'returncode': returncode,
               'outputs': [provenance.file_record(directory / name)
                           for name in filenames + ['analysis_metadata.json']]}
    (directory / 'run_manifest.json').write_text(json.dumps(receipt), encoding='utf-8')


def assert_verified(result, basis):
    assert result['errors'] == []
    assert result['validation_status'] == 'passed'
    assert result['provenance_status'] == 'verified'
    assert result['provenance_basis'] == basis
    assert result['full_validation_passed'] is True
    assert 'not scientific certification' in result['validation_scope']
    assert result['patient_multiomics_fusion_confirmed'] is False


def assert_incomplete(result):
    assert result['errors'] == []
    assert result['validation_status'] == 'incomplete_provenance'
    assert result['provenance_status'] == 'incomplete'
    assert result['provenance_basis'] is None
    assert result['full_validation_passed'] is False


def assert_invalid(result):
    assert result['errors']
    assert result['validation_status'] == 'failed'
    assert result['provenance_status'] == 'invalid'
    assert result['provenance_basis'] is None
    assert result['full_validation_passed'] is False


def test_st_native_result_hashes_need_no_run_or_metadata_self_hash(tmp_path):
    metadata, filenames = prepare(tmp_path)
    add_metadata_hashes(tmp_path, metadata, filenames)
    assert not (tmp_path / 'run_manifest.json').exists()
    assert_verified(validate_output(tmp_path), 'result_hash_inventory')


def test_gse127530_two_tables_without_result_list_use_completed_run(tmp_path):
    metadata, filenames = prepare(tmp_path, 'GSE127530')
    assert 'result_files' not in metadata
    add_run_hashes(tmp_path, filenames)
    assert_verified(validate_output(tmp_path), 'completed_run_hash_inventory')


def test_filename_only_inventory_is_explicitly_incomplete(tmp_path):
    metadata, filenames = prepare(tmp_path, 'GSE127530')
    metadata['result_files'] = filenames
    save_metadata(tmp_path, metadata)
    assert_incomplete(validate_output(tmp_path))


def test_filename_inventory_plus_completed_run_is_verified(tmp_path):
    metadata, filenames = prepare(tmp_path, 'GSE127530')
    metadata['result_files'] = filenames
    save_metadata(tmp_path, metadata)
    add_run_hashes(tmp_path, filenames)
    assert_verified(validate_output(tmp_path), 'completed_run_hash_inventory')


@pytest.mark.parametrize('version', [False, True])
def test_absent_both_receipts_preserve_legacy_structure_but_not_full_success(tmp_path, version):
    prepare(tmp_path, version=version)
    assert_incomplete(validate_output(tmp_path))


@pytest.mark.parametrize('records', [None, [], 'not-a-list'])
def test_malformed_present_metadata_inventory_is_invalid(tmp_path, records):
    metadata, _ = prepare(tmp_path)
    metadata['result_files'] = records
    save_metadata(tmp_path, metadata)
    assert_invalid(validate_output(tmp_path))


@pytest.mark.parametrize('receipt', [None, '{', {'status': 'completed', 'returncode': 0, 'outputs': None},
                                    {'status': 'completed', 'returncode': 0, 'outputs': []}])
def test_malformed_present_run_is_invalid(tmp_path, receipt):
    prepare(tmp_path)
    payload = receipt if isinstance(receipt, str) else json.dumps(receipt)
    (tmp_path / 'run_manifest.json').write_text(payload, encoding='utf-8')
    assert_invalid(validate_output(tmp_path))


def test_valid_metadata_inventory_does_not_override_invalid_run(tmp_path):
    metadata, filenames = prepare(tmp_path)
    add_metadata_hashes(tmp_path, metadata, filenames)
    add_run_hashes(tmp_path, filenames, status='failed', returncode=1)
    assert_invalid(validate_output(tmp_path))


def test_valid_run_does_not_override_invalid_metadata_inventory(tmp_path):
    metadata, filenames = prepare(tmp_path)
    metadata['result_files'] = None
    save_metadata(tmp_path, metadata)
    add_run_hashes(tmp_path, filenames)
    assert_invalid(validate_output(tmp_path))


def test_valid_hashes_do_not_override_invalid_table(tmp_path):
    metadata, filenames = prepare(tmp_path)
    table = pd.read_csv(tmp_path / filenames[0], sep='\t')
    table.loc[0, 'standard_error'] = -0.1
    table.to_csv(tmp_path / filenames[0], sep='\t', index=False, compression='gzip')
    add_metadata_hashes(tmp_path, metadata, filenames)
    result = validate_output(tmp_path)
    assert result['errors']
    assert result['provenance_status'] == 'verified'
    assert result['validation_status'] == 'failed'
    assert result['full_validation_passed'] is False


@pytest.mark.parametrize(('case', 'expected_exit'), [('st_valid', 0), ('gse_valid', 0),
                                                    ('absent', 2), ('filename_only', 2),
                                                    ('invalid_inventory', 1)])
def test_actual_cli_validate_exit_gate(tmp_path, monkeypatch, capsys, case, expected_exit):
    metadata, filenames = prepare(tmp_path, 'GSE127530' if case == 'gse_valid' else 'ST001257')
    if case == 'st_valid':
        add_metadata_hashes(tmp_path, metadata, filenames)
    elif case == 'gse_valid':
        add_run_hashes(tmp_path, filenames)
    elif case == 'filename_only':
        metadata['result_files'] = filenames
        save_metadata(tmp_path, metadata)
    elif case == 'invalid_inventory':
        metadata['result_files'] = None
        save_metadata(tmp_path, metadata)

    # Inject only the exact tested validator into main's ordinary relative import.
    # No full SciPy-dependent package/module import is claimed by this harness.
    package = types.ModuleType('nutriomics_atlas')
    package.__path__ = []
    validation_module = types.ModuleType('nutriomics_atlas.validation')
    validation_module.validate_output = validate_output
    monkeypatch.setitem(sys.modules, 'nutriomics_atlas', package)
    monkeypatch.setitem(sys.modules, 'nutriomics_atlas.validation', validation_module)
    try:
        cli_main(['validate', '--input-dir', str(tmp_path)])
        actual_exit = 0
    except SystemExit as exc:
        actual_exit = exc.code
    assert actual_exit == expected_exit
    report = json.loads(capsys.readouterr().out)
    assert report['full_validation_passed'] is (expected_exit == 0)
