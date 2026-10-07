"""Strict native BindingDB parsing, with rejection counts and raw provenance."""
import csv
import hashlib
import io
import math
import re
import zipfile
from collections import Counter
from functools import lru_cache
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from .io import sha256, utc_now, write_json

RDLogger.DisableLog("rdApp.warning")
AA = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ")
NUMBER = re.compile(r"^[+]?\d*\.?\d+(?:[eE][+-]?\d+)?$")
UNCERTAIN_CONSTRUCT = re.compile(r"\b[A-Z]\d{2,}[A-Z]\b|\bmut(?:ant|ation)\b|truncat|fusion|\[\d+\s*[-:]\s*\d+\]",re.IGNORECASE)


def uncertain_construct(name):
    return bool(UNCERTAIN_CONSTRUCT.search(name))


def to_pkd(value, endpoint="Kd", unit="nM", relation="="):
    if endpoint != "Kd":
        raise ValueError("Only Kd is admitted to this regression task")
    text = str(value).strip()
    if text.startswith("="):
        text = text[1:].strip()
    if relation != "=" or not NUMBER.fullmatch(text):
        raise ValueError("Censored or nonnumeric affinity")
    number = float(text)
    multipliers = {"nM": 1e-9, "uM": 1e-6, "µM": 1e-6, "M": 1.0, "pM": 1e-12}
    if unit not in multipliers or not math.isfinite(number) or number <= 0:
        raise ValueError("Invalid unit or nonpositive affinity")
    return -math.log10(number * multipliers[unit])


@lru_cache(maxsize=300000)
def structure(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError("Invalid SMILES")
    # Retain full stereo and disconnected components. No silent salt/tautomer merging.
    canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=molecule, includeChirality=False)
    if not scaffold:
        scaffold = "acyclic:" + canonical
    return canonical, Chem.MolToInchiKey(molecule), scaffold


def rows(path):
    csv.field_size_limit(64 * 1024 * 1024)
    if str(path).endswith(".zip"):
        archive = zipfile.ZipFile(path)
        members = [name for name in archive.namelist() if name.endswith((".tsv", ".txt"))]
        if len(members) != 1:
            raise ValueError(f"Expected one TSV archive member, found {members}")
        with archive, archive.open(members[0]) as binary, io.TextIOWrapper(binary, encoding="utf-8", errors="replace", newline="") as stream:
            yield from csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)
    else:
        with open(path, encoding="utf-8", errors="replace", newline="") as stream:
            yield from csv.reader(stream, delimiter="\t", quoting=csv.QUOTE_NONE)


def assay_lookup(path):
    if path is None:
        return {}
    lookup = {}
    for row in rows(path):
        if len(row) >= 2 and row[0].strip().isdigit():
            lookup.setdefault(row[0].strip(), set()).add(row[1].strip())
    return {key: ";".join(sorted(value)) for key, value in lookup.items()}


def prepare(raw, out, assay_map=None, species=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    mapping = assay_lookup(assay_map)
    iterator = rows(raw)
    header = next(iterator)
    header[0] = header[0].lstrip("\ufeff")
    indexes = {}
    for index, name in enumerate(header):
        indexes.setdefault(name, index)
    # Native 202609 schema numbers chain columns, while older exports do not.
    for name in list(indexes):
        if name.endswith(" 1") and ("Target Chain" in name or "Target Chain Sequence" in name):
            indexes.setdefault(name[:-2],indexes[name])
    needed = ["Kd (nM)", "Ligand SMILES", "BindingDB Target Chain Sequence"]
    for column in needed:
        if column not in indexes:
            raise ValueError(f"Native schema missing {column}; columns: {header}")
    sequence_index = indexes["BindingDB Target Chain Sequence"]
    counts = Counter()
    kept = []
    for raw_row in iterator:
        counts["raw_rows"] += 1
        def get(name):
            index = indexes.get(name)
            return raw_row[index].strip() if index is not None and index < len(raw_row) else ""
        kd_raw = get("Kd (nM)")
        if not kd_raw:
            counts["no_kd"] += 1
            continue
        try:
            pkd = to_pkd(kd_raw)
        except ValueError:
            counts["censored_or_invalid_kd"] += 1
            continue
        chains = get("Number of Protein Chains in Target (>1 implies a multichain complex)")
        if chains and chains != "1":
            counts["multichain_target"] += 1
            continue
        if get("BindingDB Target Chain Sequence 2"):
            counts["multichain_target"] += 1
            continue
        target_name = get("Target Name")
        if uncertain_construct(target_name):
            counts["explicit_mutant_truncation_or_fusion"] += 1
            continue
        sequence = raw_row[sequence_index].strip().upper() if sequence_index < len(raw_row) else ""
        if len(sequence) < 20 or not set(sequence).issubset(AA):
            counts["invalid_sequence"] += 1
            continue
        organism = get("Target Source Organism According to Curator or DataSource")
        if species and organism.casefold() != species.casefold():
            counts["other_species"] += 1
            continue
        try:
            smiles, compound_id, scaffold = structure(get("Ligand SMILES"))
        except ValueError:
            counts["invalid_smiles"] += 1
            continue
        rsid = get("BindingDB Reactant_set_id")
        doi, pmid = get("Article DOI").lower(), get("PMID")
        patent, aid = get("Patent Number"), get("PubChem AID")
        assay_id = mapping.get(rsid, "")
        # All identifiers retained, with a conservative paper-level primary grouping.
        if doi:
            source_group = "doi:" + doi
        elif pmid:
            source_group = "pmid:" + pmid
        elif patent:
            source_group = "patent:" + patent
        elif aid:
            source_group = "pubchem_aid:" + aid
        elif assay_id:
            source_group = "assay:" + assay_id
        else:
            counts["missing_independent_source_group"] += 1
            continue
        target_id = hashlib.sha256(sequence.encode()).hexdigest()[:24]
        pair_id = compound_id + ":" + target_id
        measurement_id = hashlib.sha256("|".join([pair_id, source_group, assay_id, "Kd", kd_raw, get("pH"), get("Temp (C)")]).encode()).hexdigest()[:24]
        record = {"measurement_id": measurement_id, "pair_id": pair_id, "compound_id": compound_id,
                  "target_id": target_id, "smiles": smiles, "protein_sequence": sequence,
                  "scaffold": scaffold, "target_name": target_name, "organism": organism,
                  "endpoint": "Kd", "relation": "=", "unit": "nM", "value_raw": kd_raw, "value_nm": float(kd_raw.lstrip("=").strip()), "pkd": pkd,
                  "source_group": source_group, "doi": doi, "pmid": pmid, "patent": patent, "pubchem_aid": aid,
                  "assay_id": assay_id, "reactant_set_id": rsid, "curation_source": get("Curation/DataSource"),
                  "bindingdb_monomer_id": get("BindingDB MonomerID"), "bindingdb_inchikey": get("Ligand InChI Key"),
                  "ph": get("pH"), "temperature_c": get("Temp (C)"),
                  "uniprot": get("UniProt (SwissProt) Primary ID of Target Chain") or get("UniProt (TrEMBL) Primary ID of Target Chain"),
                  "publication_date": get("Publication Date") or get("Article Publication Date"),
                  "curation_date": get("Curation Date") or get("BindingDB Curation Date")}
        kept.append(record)
        if counts["raw_rows"] % 500000 == 0:
            print(dict(counts), flush=True)
    if not kept:
        raise ValueError(f"No eligible measured Kd rows: {dict(counts)}")
    frame = pd.DataFrame(kept)
    # Alias papers across missing DOI/PMID and imports, and overlapping native assays.
    from .splits import UnionFind
    aliases = []
    for record in kept:
        ids = [record["source_group"]]
        ids += [prefix+record[field] for field,prefix in (("doi","doi:"),("pmid","pmid:"),("patent","patent:"),("pubchem_aid","pubchem_aid:")) if record[field]]
        ids += ["assay:"+value for value in record["assay_id"].split(";") if value]
        aliases.append(sorted(set(ids)))
    union = UnionFind(set(value for identifiers in aliases for value in identifiers))
    for identifiers in aliases:
        for value in identifiers[1:]:
            union.union(identifiers[0],value)
    frame["source_group"] = [union.find(identifiers[0]) for identifiers in aliases]
    frame["measurement_id"] = [hashlib.sha256("|".join([row.pair_id,row.source_group,row.assay_id,"Kd",format(row.value_nm,".17g"),row.ph,row.temperature_c,row.target_name,row.organism]).encode()).hexdigest()[:24] for row in frame.itertuples()]
    before = len(frame)
    frame = frame.drop_duplicates("measurement_id").sort_values("measurement_id").reset_index(drop=True)
    counts["duplicate_measurements"] = before - len(frame)
    counts["retained_measurements"] = len(frame)
    frame.to_csv(out / "measurements.csv.gz", index=False)
    audit = {"prepared_utc": utc_now(), "raw_file": str(Path(raw).name), "raw_sha256": sha256(raw),
             "assay_map_sha256": sha256(assay_map) if assay_map else None, "counts": dict(counts),
             "compounds": frame.compound_id.nunique(), "targets": frame.target_id.nunique(),
             "source_groups": frame.source_group.nunique(), "species_filter": species,
             "endpoint": "exact positive Kd only; native nM; pKd=-log10(Kd[M]); explicitly named mutant/truncated/fusion constructs excluded",
             "columns": header, "external_validation": "not yet established; internal source-aware holdouts only",
             "limitations": ["Native UniProt sequences may omit assayed mutations/truncations/tags; protein construct remains uncertain.",
                             "Construct-name filter is conservative pattern matching, not complete experimental sequence verification; short mutation tokens may remain.",
                             "Not a food-compound cohort: food identity annotation is required before food-specific claims.",
                             "No PubChem/ChEMBL independent external assay cohort has been acquired."]}
    write_json(out / "prepare_audit.json", audit)
    return audit
