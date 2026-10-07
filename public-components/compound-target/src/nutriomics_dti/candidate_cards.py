"""Source-grounded food-analyte research cards, without a benefit classifier.

Diseases are an explicit exact-ID protocol. CTD associations retain their
native evidence kinds; CAS/CID mappings without a full key remain uncertain.
"""
from collections import Counter
import json
import re
import pandas as pd
from rdkit import Chem


DISEASES = {
    "nafld": {"disease_id": "MESH:D065626", "expected_name": "Non-alcoholic Fatty Liver Disease", "requested": 17},
    "immune_related_inflammation": {"disease_id": "MESH:D007249", "expected_name": "Inflammation", "requested": 6},
}


def tokens(value, separator="|"):
    return set(str(value).split(separator)) - {"", "nan", "None"}


def chemical_id(value):
    text = str(value).removeprefix("MESH:")
    if not re.fullmatch(r"[CD]\d+", text):
        raise ValueError("Chemical identifier outside explicit CTD MeSH namespace")
    return "MESH:" + text


def cas_checksum(value):
    if not re.fullmatch(r"[1-9]\d{1,6}-\d{2}-\d", value):
        return False
    digits = value.replace("-", "")
    return sum(int(digit) * index for index, digit in enumerate(digits[-2::-1], 1)) % 10 == int(digits[-1])


def disease_protocol(registry):
    resolved = {}
    for task, definition in DISEASES.items():
        rows = registry[registry.DiseaseID.eq(definition["disease_id"])]
        if len(rows) != 1 or rows.iloc[0].DiseaseName != definition["expected_name"]:
            raise ValueError("Exact disease ID/name must agree with frozen CTD registry")
        resolved[task] = {**definition, "registry": rows.iloc[0].to_dict(),
                          "descendants_included": False, "name_search_used": False}
    return resolved


def map_food_chemicals(panel, registry, synonyms):
    required = {"compound_id", "smiles", "pubchem_cid", "food_verified"}
    if not required <= set(panel) or panel.compound_id.duplicated().any():
        raise ValueError("Source-verified unique structural food panel required")
    if not panel.food_verified.astype(str).str.lower().eq("true").all():
        raise ValueError("Food membership is unverified")
    cid_index, key_index, cas_index = {}, {}, {}
    for row in registry.to_dict("records"):
        for key in tokens(row.get("InChIKey", "")):
            key_index.setdefault(key, []).append(row)
        for cid in tokens(row.get("PubChemCID", "")):
            cid_index.setdefault(cid.removeprefix("CID:"), []).append(row)
        for cas in tokens(row.get("CasRN", "")):
            if cas_checksum(cas):
                cas_index.setdefault(cas, []).append(row)
    maps, rejected = [], []
    for food in panel.to_dict("records"):
        key = food["compound_id"]
        mol = Chem.MolFromSmiles(food["smiles"])
        if mol is None or Chem.MolToInchiKey(mol) != key:
            raise ValueError("Food molecular identity fails full structure guard")
        cid = str(int(food["pubchem_cid"]))
        exact = key_index.get(key, [])
        cas = sorted(value for value in synonyms.get(cid, []) if cas_checksum(value))
        candidates = exact
        basis = "exact_full_inchikey"
        if not candidates:
            candidates = [row for row in cid_index.get(cid, []) if not tokens(row.get("InChIKey", ""))]
            basis = "exact_pubchem_cid_CTD_fullkey_absent"
        if not candidates:
            candidates = [row for value in cas for row in cas_index.get(value, [])
                          if not tokens(row.get("InChIKey", ""))]
            basis = "unique_source_CAS_CTD_fullkey_absent"
        unique = {chemical_id(row["ChemicalID"]): row for row in candidates}
        if len(unique) != 1:
            rejected.append({"compound_id": key, "name": food["name"], "reason": "unmapped_or_multiple_CTD_concepts",
                             "candidate_chemical_ids": sorted(unique), "pubchem_source_CAS": cas})
            continue
        chemical = next(iter(unique.values()))
        maps.append({"compound_id": key, "name": food["name"], "pubchem_cid": int(cid),
                     "ctd_chemical_id": chemical_id(chemical["ChemicalID"]), "ctd_name": chemical["ChemicalName"],
                     "ctd_cas": chemical.get("CasRN", ""), "ctd_inchikey": chemical.get("InChIKey", ""),
                     "pubchem_source_CAS": cas, "identity_basis": basis,
                     "full_stereochemical_identity_verified": basis == "exact_full_inchikey",
                     "identity_uncertainty": None if basis == "exact_full_inchikey" else
                     "CTD parent/standard CAS or CID concept does not establish native food stereochemistry, charge, salt or conjugate",
                     "food_basis": food.get("food_basis", ""), "food_membership_source": food.get("food_membership_source", ""),
                     "positive_food_records": int(food.get("positive_food_records", 0)),
                     "pubchem_url": food.get("pubchem_url", ""), "smiles": food["smiles"]})
    return maps, rejected


def build_cards(panel, chemicals, diseases, associations, synonyms, provenance):
    protocol = disease_protocol(diseases)
    mappings, rejected = map_food_chemicals(panel, chemicals, synonyms)
    associations = associations.copy()
    associations["canonical_chemical_id"] = associations.ChemicalID.map(chemical_id)
    cards, summaries = [], {}
    for task, definition in protocol.items():
        disease_id = definition["disease_id"]
        candidates = []
        for mapping in mappings:
            subset = associations[associations.canonical_chemical_id.eq(mapping["ctd_chemical_id"])
                                  & associations.DiseaseID.eq(disease_id)]
            evidence, therapeutic, all_pmids = [], set(), set()
            for row in subset.to_dict("records"):
                kinds = sorted(tokens(row["DirectEvidence"]))
                pmids = sorted(tokens(row["PubMedIDs"]))
                if not kinds or not pmids:
                    continue
                if any(not re.fullmatch(r"[1-9]\d*", pmid) for pmid in pmids):
                    raise ValueError("Invalid native PubMed identifier")
                if not set(kinds) <= {"therapeutic", "marker/mechanism"}:
                    raise ValueError("Unexpected CTD curated evidence kind")
                all_pmids.update(pmids)
                if "therapeutic" in kinds:
                    therapeutic.update(pmids)
                evidence.append({"disease_id": disease_id, "chemical_id": row["ChemicalID"],
                                 "evidence_kinds": kinds, "pmids": pmids,
                                 "pubmed_urls": [f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" for pmid in pmids],
                                 "species": "not reported in this CTD chemical-disease download",
                                 "direction": "not separately encoded in the downloaded curated association",
                                 "source": "CTD_curated_chemicals_diseases", "source_row": row})
            if not evidence:
                continue
            candidates.append({"task": task, "disease_id": disease_id,
                               "disease_name": definition["expected_name"], **mapping,
                               "therapeutic_evidence_pmids": sorted(therapeutic), "all_evidence_pmids": sorted(all_pmids),
                               "unique_therapeutic_pmids": len(therapeutic), "unique_total_pmids": len(all_pmids),
                               "evidence": evidence, "status": "source-grounded research priority; not an experimental hit or clinical benefit",
                               "efficacy_verified": False, "fresh_DTI_score_used": False,
                               "uncertainty": ["No dose, absorption or individual health outcome inferred",
                                               "Species and direction require source-paper review",
                                               "Marker/mechanism association alone is not evidence of benefit",
                                               "Analytical aglycone-equivalent food membership differs from native free-ligand exposure"],
                               "provenance": provenance})
        candidates.sort(key=lambda row: (-row["unique_therapeutic_pmids"], -row["unique_total_pmids"],
                                         not row["full_stereochemical_identity_verified"], row["compound_id"]))
        for index, card in enumerate(candidates, 1):
            card["rank"] = index
            card["requested_top_panel"] = index <= definition["requested"]
            card["ranking_rule"] = "Descending unique therapeutic PMID count, total PMID count, full-key certainty, then full InChIKey; unvalidated prioritization"
        cards.extend(candidates)
        available = len(candidates)
        summaries[task] = {"requested": definition["requested"], "supported_candidates": available,
                           "selected": min(available, definition["requested"]), "shortfall": max(0, definition["requested"] - available),
                           "numerical_candidate_count_met": available >= definition["requested"],
                           "experimental_hit_count": 0, "clinical_benefit_claim": False,
                           "therapeutic_evidence_candidates": sum(bool(card["unique_therapeutic_pmids"]) for card in candidates)}
    return {"protocol": protocol, "summaries": summaries, "cards": cards,
            "mappings": mappings, "excluded_mappings": rejected,
            "species_information": "Not encoded by native chemical-disease association; no human efficacy assumed",
            "provenance": provenance}
