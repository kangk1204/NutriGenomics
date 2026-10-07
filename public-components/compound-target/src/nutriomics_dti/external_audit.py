"""Source/pair exclusion for externally acquired affinity records.

Database names are not independence evidence. Complete source and measurement
ledgers are required before reporting an independent regression estimate.
"""
from __future__ import annotations
from collections import Counter
import math
import re
import pandas as pd
from .data import to_pkd


def aliases(row):
    output = set()
    for field, prefix in (("doi", "doi"), ("pmid", "pmid"), ("patent", "patent"),
                          ("pubchem_aid", "pubchem_aid"), ("assay_id", "assay"),
                          ("assay_chembl_id", "chembl_assay"), ("document_chembl_id", "chembl_document")):
        value = str(row.get(field, "") or "").strip()
        if value.lower() in {"", "nan", "none", "null"}:
            continue
        if field == "doi":
            value = re.sub(r"^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)", "", value, flags=re.I).lower()
        if field == "pmid":
            value = value.removesuffix(".0")
            if not value.isdigit():
                raise ValueError("Invalid PMID in source ledger")
        output.add(prefix + ":" + value)
    group = str(row.get("source_group", "") or "").strip()
    if group and group.lower() not in {"nan", "none", "null"}:
        output.add(group)
    return output


def audit(reference, external, complete_reference=False):
    required = {"compound_id", "target_id"}
    if not required <= set(reference) or not required <= set(external):
        raise ValueError("Exact compound/target identifiers are required")
    seen_pairs = set(zip(reference.compound_id, reference.target_id))
    source_ids = set()
    reference_without_source = 0
    for row in reference.to_dict("records"):
        origins = aliases(row)
        source_ids |= origins
        reference_without_source += not bool(origins)
    reasons = Counter()
    eligible = []
    decisions = []
    for index, row in enumerate(external.to_dict("records")):
        why = []
        try:
            pkd = to_pkd(row.get("value_nm"), endpoint=row.get("endpoint"), unit="nM", relation=row.get("relation"))
            if not math.isfinite(pkd):
                raise ValueError("Nonfinite pKd")
        except (TypeError, ValueError):
            why.append("unsupported_or_invalid_exact_positive_kd")
        if not row["compound_id"] or not row["target_id"]:
            why.append("missing_structure_or_sequence_identity")
        if (row["compound_id"], row["target_id"]) in seen_pairs:
            why.append("compound_target_pair_seen_in_reference")
        origins = aliases(row)
        if not origins:
            why.append("external_source_unresolved")
        if origins & source_ids:
            why.append("paper_patent_assay_source_overlap")
        reasons.update(why)
        if not why:
            eligible.append(index)
        decisions.append({"row_index": index, "activity_id": row.get("activity_id"), "reasons": why})
    # An omitted source ledger cannot be interpreted as an empty one.
    independent = bool(complete_reference and not reference_without_source and len(eligible))
    return external.iloc[eligible].copy(), {
        "reference_records": len(reference), "external_records": len(external),
        "reference_without_source_alias": reference_without_source,
        "complete_reference_ledger_verified": bool(complete_reference),
        "source_and_pair_nonoverlapping_candidates": len(eligible),
        "rejection_reason_counts": dict(reasons), "decisions": decisions,
        "independent_evaluation_eligible": independent,
        "independence_note": "Full-key/sequence/paper/assay exclusion only; similar structures, same laboratories and upstream normalization remain separate applicability checks"}
