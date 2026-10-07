"""Streaming CTD parser. CSV quoting and gzip integrity are validated through EOF."""
from __future__ import annotations
import csv
import gzip
import itertools
import json
import os
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from .sources import atomic_json, file_record, now, registry_hook, sha256

DEFAULT_FILES = ('CTD_chem_gene_ixns.csv.gz', 'CTD_curated_chemicals_diseases.csv.gz',
                 'CTD_curated_genes_diseases.csv.gz', 'CTD_chemicals.csv.gz',
                 'CTD_diseases.csv.gz', 'CTD_genes_pathways.csv.gz', 'CTD_pathways.csv.gz')
CTD_RIGHTS = 'CTD source terms: citation, source hyperlinks, use notification and periodic access; not a redistribution permission determination'
csv.field_size_limit(2 ** 27)
REQUIRED_FIELDS = {
    'CTD_chem_gene_ixns.csv.gz': {'ChemicalName','ChemicalID','GeneSymbol','GeneID','GeneForms','Organism','OrganismID','Interaction','InteractionActions','PubMedIDs'},
    'CTD_curated_chemicals_diseases.csv.gz': {'ChemicalName','ChemicalID','DiseaseName','DiseaseID','DirectEvidence','PubMedIDs'},
    'CTD_curated_genes_diseases.csv.gz': {'GeneSymbol','GeneID','DiseaseName','DiseaseID','DirectEvidence','OmimIDs','PubMedIDs'},
    'CTD_chemicals_diseases.csv.gz': {'ChemicalName','ChemicalID','DiseaseName','DiseaseID','DirectEvidence','InferenceGeneSymbol','InferenceScore','PubMedIDs'},
    'CTD_chemicals.csv.gz': {'ChemicalName','ChemicalID','PubChemCID','InChIKey'},
    'CTD_genes.csv.gz': {'GeneSymbol','GeneID','UniProtIDs'},
    'CTD_genes_pathways.csv.gz': {'GeneSymbol','GeneID','PathwayName','PathwayID'},
}


def mesh_id(value: str) -> str:
    value = value.strip()
    if re.fullmatch(r'[CD]\d+', value):
        return 'MESH:' + value
    if re.fullmatch(r'(MESH:[CD]\d+|OMIM:\d+)', value):
        return value
    raise ValueError(f'unsupported disease/chemical ID: {value!r}')


def gene_id(value: str) -> str:
    if not value.isdigit():
        raise ValueError(f'invalid NCBI Gene identifier: {value!r}')
    return 'NCBIGene:' + value


def ctd_header(handle) -> tuple[dict, str]:
    comments, fields, report_date = [], None, None
    for line in handle:
        if line.startswith('#'):
            comments.append(line.rstrip())
            if line.startswith('# Report created:'):
                report_date = line.split(':', 1)[1].strip()
            if line.startswith('# Fields:'):
                schema_line = next(handle, '')
                if not schema_line.startswith('# '):
                    raise ValueError('missing commented CTD schema')
                fields = next(csv.reader([schema_line[2:].strip()]))
                # Some official exports have a trailing empty column, retained as _trailing.
                fields = [v if v else '_trailing' for v in fields]
                comments.append(schema_line.rstrip())
        elif line.strip():
            if not fields or not report_date:
                raise ValueError('CTD report date and Fields header are required')
            return {'fields': fields, 'report_created': report_date, 'comments': comments}, line
    raise ValueError('empty CTD file or missing data rows')


def read_ctd(path: Path):
    """Yield header once, then dict records. Exhaustion validates gzip CRC/footer."""
    if path.name.endswith('.crdownload') or path.suffix not in {'.gz', '.csv'}:
        raise ValueError('only completed csv/csv.gz snapshots are supported')
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        metadata, first = ctd_header(handle)
        if len(set(metadata['fields'])) != len(metadata['fields']):
            raise ValueError('duplicate CTD schema column')
        if not REQUIRED_FIELDS.get(path.name,set()).issubset(metadata['fields']):
            raise ValueError('CTD source schema missing required fields')
        yield metadata
        for index, row in enumerate(csv.reader(itertools.chain([first], handle), strict=True), 1):
            if not row:
                continue
            if len(row) != len(metadata['fields']):
                raise ValueError(f'{path.name}: data row {index} has {len(row)} fields, expected {len(metadata["fields"])}')
            yield dict(zip(metadata['fields'], row, strict=True))


def ingest_ctd(source_dir: Path, output_dir: Path, files=DEFAULT_FILES, registry: Path | None = None, batch_size: int = 25000) -> dict:
    import pyarrow as pa
    import pyarrow.parquet as pq
    import duckdb
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / 'manifest.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {'schema_version': 1, 'source': 'CTD', 'files': {}}
    for name in files:
        if Path(name).name != name or not name.startswith('CTD_') or name.endswith('.crdownload'):
            raise ValueError('unsafe or unsupported source filename')
        source = source_dir / name
        iterator = read_ctd(source)
        header = next(iterator)
        report_date = header['report_created']
        # Official timezone abbreviations are retained; date is not filesystem mtime.
        parsed = datetime.strptime(report_date.rsplit(' ', 2)[0] + ' ' + report_date.rsplit(' ', 1)[1], '%a %b %d %H:%M:%S %Y')
        record = file_record(source, 'CTD', f'https://ctdbase.org/reports/{name}', parsed.date().isoformat(), CTD_RIGHTS)
        target = output_dir / (name.removesuffix('.gz').removesuffix('.csv') + '.parquet')
        prior = manifest['files'].get(name)
        if prior and prior['sha256'] == record['sha256'] and target.exists() and sha256(target) == prior['parquet_sha256']:
            iterator.close()
            continue
        if prior:
            raise ValueError(f'new source release requires a new output directory: {name}')
        temporary = target.with_suffix('.parquet.part')
        schema = pa.schema([(field, pa.string()) for field in header['fields']])
        rows, batch, taxons, evidence, writer = 0, [], Counter(), Counter(), None
        try:
            writer = pq.ParquetWriter(temporary, schema, compression='zstd')
            for row in iterator:
                batch.append(row); rows += 1
                taxons.update([row.get('OrganismID', row.get('organismid', ''))])
                evidence.update([row.get('DirectEvidence', '')])
                if len(batch) >= batch_size:
                    writer.write_table(pa.Table.from_pylist(batch, schema=schema)); batch.clear()
            if batch:
                writer.write_table(pa.Table.from_pylist(batch, schema=schema))
            writer.close(); writer = None
            if sha256(source) != record['sha256']:
                raise ValueError('source snapshot changed during ingestion: ' + name)
            os.replace(temporary, target)
        except BaseException:
            if writer is not None:
                writer.close()
            temporary.unlink(missing_ok=True)
            raise
        record.update({'rows': rows, 'fields': header['fields'], 'report_created': report_date,
                       'gzip_crc': 'passed' if source.suffix == '.gz' else 'not_applicable',
                       'parquet': target.name, 'parquet_sha256': sha256(target),
                       'taxon_counts': dict(taxons), 'direct_evidence_counts': dict(evidence)})
        manifest['files'][name] = record
        manifest['updated_at'] = now()
        atomic_json(manifest_path, manifest)
        registry_hook(registry, record)
    # DuckDB views remain relocatable because the catalog is rebuilt from the manifest.
    catalog = duckdb.connect(str(output_dir / 'ctd.duckdb'))
    try:
        for record in manifest['files'].values():
            name = Path(record['parquet']).stem.replace('"', '""')
            path = str((output_dir / record['parquet']).resolve()).replace("'", "''")
            catalog.execute(f'CREATE OR REPLACE VIEW "{name}" AS SELECT * FROM read_parquet(\'{path}\')')
    finally:
        catalog.close()
    return manifest
