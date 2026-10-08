from __future__ import annotations
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from .metabolomics import bh
from .provenance import sha256


def validate_output(directory: Path) -> dict:
    """Validate recorded tables and receipts without changing analysis results."""
    directory = Path(directory)
    errors = []
    tables = {}
    try:
        metadata = json.loads((directory / 'analysis_metadata.json').read_text(encoding='utf-8'))
        if not isinstance(metadata, dict):
            raise ValueError('metadata must be an object')
    except (OSError, ValueError) as exc:
        return {'accession': None, 'independent_people': None, 'tables': {},
                'errors': ['invalid analysis metadata: ' + str(exc)],
                'patient_multiomics_fusion_confirmed': False}

    accession = metadata.get('accession')
    cohort = metadata.get('independent_people')
    try:
        cohort_number = float(cohort)
        cohort_valid = not isinstance(cohort, bool) and np.isfinite(cohort_number) and cohort_number > 0 and cohort_number.is_integer()
    except (TypeError, ValueError):
        cohort_valid = False
    if not cohort_valid:
        errors.append('invalid independent cohort count')
    if not isinstance(accession, str) or not accession.strip():
        errors.append('invalid accession')
    expected = ['metabolite_effects.tsv.gz'] if accession == 'ST001257' else ['gene_effects.tsv.gz', 'pathway_effects.tsv.gz']

    for name in expected:
        path = directory / name
        if not path.is_file():
            errors.append('missing result ' + name)
            continue
        try:
            table = pd.read_csv(path, sep='\t')
        except (OSError, ValueError, pd.errors.ParserError) as exc:
            errors.append('unreadable result ' + name + ': ' + str(exc))
            continue
        required = {'feature_id', 'contrast', 'effect_estimate', 'p_value', 'q_value', 'n_people'}
        if name != 'metabolite_effects.tsv.gz':
            required.add('n_samples')
        if not required.issubset(table):
            errors.append('missing result fields ' + name)
            continue
        if table.empty:
            errors.append('empty result ' + name)
            continue
        for column in ('feature_id', 'contrast'):
            if table[column].isna().any() or table[column].astype(str).str.strip().eq('').any():
                errors.append('invalid result identity ' + column + ' ' + name)
        if table.duplicated(['feature_id', 'contrast']).any():
            errors.append('duplicate feature-contrast ' + name)

        numeric = {}
        numeric_invalid = set()
        columns = {'effect_estimate', 'p_value', 'q_value', 'n_people', 'n_samples',
                   'q_value_study_family', 'standard_error', 'confidence_low', 'confidence_high'}
        for column in sorted(columns.intersection(table.columns)):
            values = pd.to_numeric(table[column], errors='coerce').to_numpy(dtype=float)
            numeric[column] = values
            if (table[column].notna().to_numpy() & np.isnan(values)).any():
                numeric_invalid.add(column)
                errors.append('nonnumeric result field ' + column + ' ' + name)
        p = numeric['p_value']
        estimable = np.isfinite(p) & (p >= 0) & (p <= 1)
        for column in ('p_value', 'q_value', 'q_value_study_family'):
            if column in numeric:
                values = numeric[column]
                if (np.isinf(values) | (np.isfinite(values) & ((values < 0) | (values > 1)))).any():
                    errors.append('invalid p/q range ' + column + ' ' + name)
        if 'q_value_study_family' not in table:
            errors.append('primary study-level testing family absent ' + name)
        elif 'p_value' not in numeric_invalid and not (np.isinf(p) | (np.isfinite(p) & ((p < 0) | (p > 1)))).any():
            family = numeric['q_value_study_family']
            if not np.allclose(family, bh(p), rtol=0, atol=1e-10, equal_nan=True):
                errors.append('BH family inconsistent ' + name)
            if not np.allclose(numeric['q_value'], family, rtol=0, atol=1e-10, equal_nan=True):
                errors.append('primary q alias inconsistent ' + name)

        status = table['status'].fillna('').astype(str) if 'status' in table else pd.Series('', index=table.index)
        nonestimable = status.isin({'not_estimable', 'insufficient_pairs_or_constant_delta'}).to_numpy()
        if (np.isnan(p) & ~nonestimable).any():
            errors.append('missing p without explicit non-estimable status ' + name)
        if (estimable & nonestimable).any():
            errors.append('non-estimable status has a p value ' + name)
        effect = numeric['effect_estimate']
        if (~np.isfinite(effect[estimable])).any():
            errors.append('estimable effect absent or nonfinite ' + name)

        people = numeric['n_people']
        if (~np.isfinite(people) | (people < 0) | (people != np.floor(people))).any() or (people[estimable] <= 0).any():
            errors.append('invalid person count ' + name)
        if cohort_valid and (people > cohort_number).any():
            errors.append('person count exceeds independent cohort ' + name)
        if 'n_samples' in numeric:
            samples = numeric['n_samples']
            if (~np.isfinite(samples) | (samples < 0) | (samples != np.floor(samples))).any() or (samples[estimable] <= 0).any() or (samples < people).any():
                errors.append('invalid sample count ' + name)

        if name != 'pathway_effects.tsv.gz':
            uncertainty = {'confidence_low', 'confidence_high', 'standard_error'}
            if not uncertainty.issubset(table):
                errors.append('gene/metabolite uncertainty absent ' + name)
            else:
                se = numeric['standard_error']
                low = numeric['confidence_low']
                high = numeric['confidence_high']
                if (np.isfinite(se) & (se < 0)).any() or (~np.isfinite(se[estimable])).any() or (se[estimable] <= 0).any():
                    errors.append('invalid standard error ' + name)
                if (~np.isfinite(low[estimable])).any() or (~np.isfinite(high[estimable])).any():
                    errors.append('estimable confidence interval absent or nonfinite ' + name)
                bounded = np.isfinite(low) & np.isfinite(high)
                if (bounded & (low > high)).any():
                    errors.append('reversed confidence interval ' + name)
                if (bounded & np.isfinite(effect) & ((effect < low) | (effect > high))).any():
                    errors.append('confidence interval excludes effect ' + name)
        tables[name] = {'rows': len(table), 'features': table.feature_id.nunique(),
                        'contrasts': table.contrast.nunique(), 'estimable_tests': int(estimable.sum()),
                        'q_below_005': int((numeric['q_value'] < .05).sum()), 'sha256': sha256(path)}

    def inspect_inventory(records, label, require_hash=False):
        names = set()
        hashed = set()
        if not isinstance(records, list):
            errors.append('invalid ' + label + ' inventory')
            return names, hashed
        for record in records:
            if isinstance(record, str) and not require_hash:
                filename = re.split(r'[\\/]', record)[-1]
                digest = None
            elif isinstance(record, dict) and isinstance(record.get('path'), str):
                filename = re.split(r'[\\/]', record['path'])[-1]
                digest = record.get('sha256')
                if not isinstance(digest, str) or not re.fullmatch(r'[0-9a-fA-F]{64}', digest):
                    errors.append('invalid ' + label + ' hash ' + filename)
                    digest = None
            else:
                errors.append('invalid ' + label + ' record')
                continue
            if not filename or filename in {'.', '..'}:
                errors.append('invalid ' + label + ' filename')
                continue
            if filename in names:
                errors.append('duplicate ' + label + ' filename ' + filename)
            names.add(filename)
            candidate = directory / filename
            if not candidate.is_file():
                errors.append('missing ' + label + ' file ' + filename)
            elif digest is not None:
                hashed.add(filename)
                if sha256(candidate) != digest.lower():
                    errors.append(label + ' provenance hash mismatch ' + filename)
        return names, hashed

    if 'result_files' in metadata:
        records = metadata['result_files']
        names, hashed = inspect_inventory(records, 'result')
        if not set(expected).issubset(names):
            errors.append('result inventory omits expected tables')
        if isinstance(records, list) and any(isinstance(record, dict) for record in records) and not set(expected).issubset(hashed):
            errors.append('result hash inventory omits expected tables')
    if 'auxiliary_files' in metadata:
        inspect_inventory(metadata['auxiliary_files'], 'auxiliary')

    run = directory / 'run_manifest.json'
    if run.exists():
        try:
            receipt = json.loads(run.read_text(encoding='utf-8'))
            if not isinstance(receipt, dict):
                raise ValueError('run manifest must be an object')
        except (OSError, ValueError) as exc:
            errors.append('invalid run manifest: ' + str(exc))
        else:
            if receipt.get('status') != 'completed' or type(receipt.get('returncode')) is not int or receipt['returncode'] != 0:
                errors.append('R run is not completed')
            _, hashed = inspect_inventory(receipt.get('outputs'), 'output', require_hash=True)
            if not (set(expected) | {'analysis_metadata.json'}).issubset(hashed):
                errors.append('output hash inventory omits required results or metadata')
    return {'accession': accession, 'independent_people': cohort, 'tables': tables,
            'errors': errors, 'patient_multiomics_fusion_confirmed': False}
