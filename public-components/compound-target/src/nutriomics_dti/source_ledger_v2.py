"""Canonical source aliases and exact native-Kd identity guards.

The ledger includes every admitted measurement, not just a food subset.
Source closure conservatively joins papers, patents and assay aliases.
"""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from .data import AA, structure, to_pkd
from .io import sha256, utc_now, write_json
from .splits import UnionFind

KEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
FIELDS = {"doi": "doi", "pmid": "pmid", "patent": "patent", "pubchem_aid": "pubchem_aid", "assay_id": "assay", "assay_chembl_id": "chembl_assay", "document_chembl_id": "chembl_document"}


def canonical_aliases(row):
    output = set()
    for field, prefix in FIELDS.items():
        value = str(row.get(field, "") or "").strip()
        if value.lower() in {"", "nan", "none", "null"}:
            continue
        for item in value.split(";"):
            item = item.strip()
            if not item:
                continue
            if field == "doi":
                item = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", item, flags=re.I).lower()
            elif field in {"pmid", "pubchem_aid"}:
                if field == "pubchem_aid":
                    item = re.sub(r"^aid\s*:?\s*", "", item, flags=re.I)
                item = item.removesuffix(".0")
                if not item.isdigit():
                    raise ValueError(f"Invalid {field}: {item}")
            elif field == "patent":
                item = re.sub(r"[\s-]", "", item).upper()
            output.add(prefix + ":" + item)
    # Existing union component is retained, but never replaces real source IDs.
    group = str(row.get("source_group", "") or "").strip()
    if group and group.lower() not in {"nan", "none", "null"}:
        output.add("component:" + group)
    return sorted(output)


def native_guards(frame):
    required = {"measurement_id", "pair_id", "compound_id", "target_id", "smiles", "protein_sequence", "endpoint", "relation", "unit", "value_nm", "pkd", "source_group"}
    if not required <= set(frame):
        raise ValueError(f"Native Kd columns missing: {sorted(required-set(frame))}")
    if frame.measurement_id.duplicated().any():
        raise ValueError("Duplicate measurement identifiers")
    for row in frame.to_dict("records"):
        if row["unit"] != "nM":
            raise ValueError("Only native nM labels admitted")
        expected = to_pkd(row["value_nm"], endpoint=row["endpoint"], relation=row["relation"])
        if not np.isfinite(float(row["pkd"])) or abs(expected-float(row["pkd"])) > 1e-8:
            raise ValueError("pKd is inconsistent with exact native positive Kd")
        _, key, _ = structure(row["smiles"])
        seq = str(row["protein_sequence"])
        tid = hashlib.sha256(seq.encode()).hexdigest()[:24]
        if not KEY.fullmatch(str(row["compound_id"])) or key != row["compound_id"]:
            raise ValueError("Full stereochemical InChIKey/SMILES mismatch")
        if len(seq) < 20 or not set(seq) <= AA or tid != row["target_id"]:
            raise ValueError("Full protein sequence/target SHA mismatch")
        if row["pair_id"] != key + ":" + tid:
            raise ValueError("Pair identifier mismatch")
        if not any(not value.startswith("component:") for value in canonical_aliases(row)):
            raise ValueError("Real paper/patent/assay source absent")


def close_aliases(frame):
    """Return a canonical union ledger; aliases cannot straddle a split."""
    frame = frame.copy()
    aliases = [canonical_aliases(row) for row in frame.to_dict("records")]
    if any(not ids for ids in aliases):
        raise ValueError("A measurement has no source alias")
    union = UnionFind({value for ids in aliases for value in ids})
    for ids in aliases:
        for value in ids[1:]:
            union.union(ids[0], value)
    frame["source_group"] = [union.find(ids[0]) for ids in aliases]
    frame["source_aliases_json"] = [json.dumps(ids, separators=(",", ":")) for ids in aliases]
    return frame


def build(prepared, out, prepare_audit=None):
    out = Path(out)
    if (out/"ledger.json").exists():
        raise FileExistsError("Frozen ledger already exists")
    frame = pd.read_csv(prepared, keep_default_na=False)
    native_guards(frame)
    frame = close_aliases(frame)
    input_records = len(frame)
    frame["input_measurement_id"] = frame.measurement_id
    # Newly discovered aliases must not turn one imported measurement into two
    # independent observations. Preserve native source and construct conditions.
    canonical_ids = []
    for row in frame.to_dict("records"):
        fields = [row["pair_id"],row["source_group"],str(row.get("assay_id","")),"Kd",format(float(row["value_nm"]),".17g")]
        fields += [str(row.get(field,"")) for field in ("ph","temperature_c","target_name","organism")]
        canonical_ids.append(hashlib.sha256("|".join(fields).encode()).hexdigest()[:24])
    frame["measurement_id"] = canonical_ids
    frame = frame.drop_duplicates("measurement_id").sort_values("measurement_id").reset_index(drop=True)
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out/"measurements.csv.gz", index=False)
    fields = ["measurement_id", "pair_id", "compound_id", "target_id", "source_group", "source_aliases_json", "doi", "pmid", "patent", "pubchem_aid", "assay_id", "reactant_set_id", "curation_source", "endpoint", "relation", "unit", "value_nm"]
    frame[[field for field in fields if field in frame]].to_csv(out/"source_ledger.csv.gz", index=False)
    prior = json.loads(Path(prepare_audit).read_text()) if prepare_audit else {}
    if prior and int(prior.get("counts", {}).get("retained_measurements", -1)) != input_records:
        raise ValueError("Prepare audit retained count disagrees with full ledger")
    full_key_equal = frame.bindingdb_inchikey.eq(frame.compound_id) if "bindingdb_inchikey" in frame else pd.Series(False,index=frame.index)
    report = {"created_utc": utc_now(), "input_sha256": sha256(prepared), "prepared_measurements_sha256": sha256(out/"measurements.csv.gz"), "ledger_sha256": sha256(out/"source_ledger.csv.gz"), "records": len(frame), "input_records":input_records, "additional_alias_duplicate_measurements":input_records-len(frame), "pairs": int(frame.pair_id.nunique()), "compounds": int(frame.compound_id.nunique()), "targets": int(frame.target_id.nunique()), "source_groups": int(frame.source_group.nunique()), "native_inchikey_full_match_records":int(full_key_equal.sum()), "native_inchikey_unresolved_or_mismatch_records":int((~full_key_equal).sum()), "raw_sha256": prior.get("raw_sha256"), "assay_map_sha256": prior.get("assay_map_sha256"), "complete_current_release_ledger": bool(prior), "historical_reference_verified": False, "endpoint": "exact positive native Kd nM", "note": "Full current release ledger does not prove complete coverage of historical training sources; model compound IDs derive from supplied isomeric SMILES. Native full-key identity uncertainty remains explicitly reported."}
    write_json(out/"ledger.json", report)
    return report


def food_external_guard(reference, external, verified_food_keys, complete_reference=False):
    """Fail closed until historical ledger and minimum food cohort are sufficient."""
    native_guards(reference)
    native_guards(external)
    keys = set(verified_food_keys)
    if any(not KEY.fullmatch(value) for value in keys):
        raise ValueError("Food panel requires full InChIKey")
    external = external[external.compound_id.isin(keys)].copy()
    prior_pairs = set(reference.pair_id)
    prior_aliases = {alias for row in reference.to_dict("records") for alias in canonical_aliases(row)}
    retained, reasons = [], []
    for index, row in external.iterrows():
        why = []
        if row.pair_id in prior_pairs:
            why.append("pair_seen")
        if set(canonical_aliases(row)) & prior_aliases:
            why.append("paper_patent_assay_seen")
        if not complete_reference:
            why.append("historical_ledger_incomplete")
        reasons.append({"measurement_id": row.measurement_id, "reasons": why})
        if not why:
            retained.append(index)
    result = external.loc[retained].copy()
    count = {"n": len(result), "compounds": int(result.compound_id.nunique()), "targets": int(result.target_id.nunique()), "source_groups": int(result.source_group.nunique())}
    minimum = {"n": 30, "compounds": 5, "targets": 3, "source_groups": 3}
    eligible = bool(complete_reference and all(count[key] >= value for key, value in minimum.items()))
    return result, {"historical_ledger_complete": bool(complete_reference), "food_records_before_exclusion": len(external), "retained": count, "minimum": minimum, "reportable": eligible, "status": "independent_food_evaluation_eligible" if eligible else "unevaluated", "decisions": reasons}
