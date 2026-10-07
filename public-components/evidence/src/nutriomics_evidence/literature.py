"""Scoped dictionary mentions and chemical-gene extraction from frozen PubMed.

Lexical mentions are candidates, not validated NER or experimental binding.
Clinical exposures/outcomes are never converted to genes or DrugProt relations.
"""
from __future__ import annotations
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from .ctd import read_ctd, CTD_RIGHTS
from .sources import atomic_json, now, sha256
from .extraction import rule_extract, sentence_spans, parse_qwen, write_relations

DISEASE_IDS = {"D006973", "D000075222", "MESH:D006973", "MESH:D000075222"}
TOPICS = {"HTN", "DASH", "omega3", "vitD", "caffeine", "cocoa_flavanol"}
RESERVED_OUTCOMES = {
    "bp", "sbp", "dbp", "htn", "dash", "blood pressure", "systolic blood pressure",
    "diastolic blood pressure", "hypertension", "arterial pressure", "dna", "dnam", "rna", "mrna",
}


def _stream(path):
    iterator = read_ctd(Path(path))
    metadata = next(iterator)
    return metadata, iterator


def ctd_dictionary(directory):
    directory = Path(directory)
    selected, files = {}, {}
    for filename, kind in (("CTD_curated_chemicals_diseases.csv.gz", "CHEMICAL"),
                           ("CTD_curated_genes_diseases.csv.gz", "GENE")):
        path = directory / filename
        metadata, iterator = _stream(path)
        count = 0
        for row in iterator:
            count += 1
            if row["DiseaseID"] not in DISEASE_IDS:
                continue
            raw = row["ChemicalID"] if kind == "CHEMICAL" else row["GeneID"]
            identifier = (raw if raw.startswith("MESH:") else "MESH:" + raw) if kind == "CHEMICAL" else "NCBIGene:" + raw
            text = row["ChemicalName"] if kind == "CHEMICAL" else row["GeneSymbol"]
            record = selected.setdefault(identifier, {"id": identifier, "type": kind, "aliases": []})
            alias = {"text": text, "case_sensitive": kind == "GENE", "basis": "CTD canonical chemical name" if kind == "CHEMICAL" else "CTD canonical gene symbol"}
            if alias not in record["aliases"]:
                record["aliases"].append(alias)
        files[filename] = {"sha256": sha256(path), "bytes": path.stat().st_size, "rows": count,
                           "header": metadata, "gzip_crc": "passed through EOF", "official_url": "https://ctdbase.org/reports/" + filename}
    gene_path = directory / "CTD_genes.csv.gz"
    if gene_path.exists():
        metadata, iterator = _stream(gene_path)
        matched, count = set(), 0
        for row in iterator:
            count += 1
            identifier = "NCBIGene:" + row["GeneID"]
            if identifier not in selected:
                continue
            matched.add(identifier)
            names = [(row.get("GeneName", ""), False, "CTD gene name")]
            names.extend((alias, not (" " in alias or alias.islower()), "CTD gene synonym")
                         for alias in row.get("Synonyms", "").split("|"))
            for text, case_sensitive, basis in names:
                if text:
                    alias = {"text": text, "case_sensitive": case_sensitive, "basis": basis}
                    if alias not in selected[identifier]["aliases"]:
                        selected[identifier]["aliases"].append(alias)
        files[gene_path.name] = {"sha256": sha256(gene_path), "bytes": gene_path.stat().st_size,
                               "rows": count, "selected_gene_ids_with_dictionary": len(matched),
                               "header": metadata, "gzip_crc": "passed through EOF",
                               "official_url": "https://ctdbase.org/reports/" + gene_path.name}
    return {"created_utc": now(), "provider": "CTD", "files": files,
            "scope_disease_ids": sorted(DISEASE_IDS), "records": sorted(selected.values(), key=lambda row: row["id"]),
            "rights": CTD_RIGHTS,
            "species_policy": "Selected NCBI gene IDs from the scoped graph; document species and clinical status are not inferred from a dictionary match"}


class MentionDictionary:
    def __init__(self, records):
        aliases = defaultdict(list)
        self.excluded = []
        for record in records:
            if record["type"] not in {"CHEMICAL", "GENE"}:
                raise ValueError("Only chemical and gene lexical candidates are permitted")
            for alias in record["aliases"]:
                text = alias["text"].strip()
                # The source synonym field contains common words and legacy
                # aliases (e.g. TH: The; LEPR: diabetes; TIMP1: EPA). These are
                # preserved in the source snapshot, but automatic gene mention
                # matching uses only canonical symbols and gene names.
                if record["type"] == "GENE" and alias["basis"] == "CTD gene synonym":
                    self.excluded.append({"id": record["id"], "alias": text,
                                          "reason": "unreviewed_source_synonym_not_automatic_gene_identity"})
                    continue
                if not text or len(text) < 3 or not re.search(r"[A-Za-z]", text):
                    self.excluded.append({"id": record["id"], "alias": text, "reason": "short_or_nonlexical_alias"})
                    continue
                if record["type"] == "GENE" and text.casefold() in RESERVED_OUTCOMES:
                    self.excluded.append({"id": record["id"], "alias": text, "reason": "clinical_outcome_or_generic_nucleic_acid_not_gene_mention"})
                    continue
                aliases[text.casefold()].append({**alias, "text": text, "id": record["id"], "type": record["type"]})
        self.aliases, self.ambiguous = {}, {}
        patterns = []
        for key, entries in aliases.items():
            ids = sorted({entry["id"] for entry in entries})
            if len(ids) > 1:
                self.ambiguous[key] = {"ids": ids, "types": sorted({entry["type"] for entry in entries})}
            self.aliases[key] = entries
            for text, case_sensitive in sorted({(entry["text"], entry["case_sensitive"]) for entry in entries}):
                expression = re.escape(text)
                patterns.append((len(text), expression if case_sensitive else "(?ai:" + expression + ")"))
        patterns.sort(key=lambda item: (-item[0], item[1]))
        self.pattern = re.compile(r"(?<!\w)(?:" + "|".join(expression for _, expression in patterns) + r")(?!\w)") if patterns else None

    def annotate(self, text):
        entities, ledger = [], []
        if self.pattern is None:
            return entities, ledger
        for match in self.pattern.finditer(text):
            token = match.group()
            entries = [entry for entry in self.aliases[token.casefold()]
                       if not entry["case_sensitive"] or entry["text"] == token]
            ids = sorted({entry["id"] for entry in entries})
            if token.casefold() == "lead" and any(entry["type"] == "CHEMICAL" for entry in entries):
                ledger.append({"start": match.start(), "end": match.end(), "text": token,
                               "candidate_ids": ids, "reason": "polysemous_chemical_name_requires_context_review"})
                continue
            if len(ids) != 1:
                ledger.append({"start": match.start(), "end": match.end(), "text": token,
                               "candidate_ids": ids, "reason": "ambiguous_dictionary_identity"})
                continue
            kinds = {entry["type"] for entry in entries}
            if len(kinds) != 1:
                raise ValueError("One dictionary ID has inconsistent lexical types")
            entities.append({"id": "T" + str(len(entities) + 1), "type": kinds.pop(), "start": match.start(),
                             "end": match.end(), "text": token, "normalized_ids": ids,
                             "identity_basis": sorted({entry["basis"] for entry in entries}),
                             "mention_status": "lexical_dictionary_candidate_requires_review"})
        return entities, ledger


def documents_from_articles(path, dictionary):
    documents, seen = [], set()
    with Path(path).open(encoding="utf-8") as stream:
        for line in stream:
            article = json.loads(line)
            pmid = article.get("pmid")
            if not isinstance(pmid, str) or not pmid.isdigit() or pmid in seen:
                raise ValueError("Duplicate/malformed PMID in frozen PubMed corpus")
            seen.add(pmid)
            if not set(article.get("topics", [])).issubset(TOPICS) or not article.get("topics"):
                raise ValueError("Article is outside the six approved query topics")
            title, abstract = article.get("title", ""), article.get("abstract", "")
            if not isinstance(title, str) or not isinstance(abstract, str):
                raise ValueError("Original title/abstract text must be literal strings")
            text = title + " " + abstract
            entities, ambiguity = dictionary.annotate(text)
            kinds = {entity["type"] for entity in entities}
            eligible = {"CHEMICAL", "GENE"}.issubset(kinds)
            same_sentence_pairs = 0
            for start, end in sentence_spans(text):
                inside = [entity for entity in entities if start <= entity["start"] and entity["end"] <= end]
                same_sentence_pairs += sum(entity["type"] == "CHEMICAL" for entity in inside) * sum(entity["type"] == "GENE" for entity in inside)
            documents.append({"pmid": pmid, "text": text, "entities": entities,
                              "ambiguity_ledger": ambiguity, "qwen_eligible": eligible,
                              "abstention_reason": None if eligible else "no_unambiguous_chemical_gene_document_comention",
                              "same_sentence_candidate_pairs": same_sentence_pairs,
                              "source": {key: article.get(key) for key in ("doi", "pmcid", "source_url", "topics", "publication_types", "mesh", "record_type")},
                              "text_policy": "literal title +one space+ literal abstract; original offsets; no normalized-text offset repair",
                              "text_sha256": hashlib.sha256(text.encode()).hexdigest()})
    return documents


def annotate_relations(relations, document):
    entities = {entity["id"]: entity for entity in document["entities"]}
    return [{**row, "subject_candidate_ids": entities[row["arg1"]]["normalized_ids"],
             "object_candidate_ids": entities[row["arg2"]]["normalized_ids"],
             "subject_text": entities[row["arg1"]]["text"], "object_text": entities[row["arg2"]]["text"],
             "source_url": document["source"]["source_url"], "source_text_sha256": document["text_sha256"],
             "mention_status": "dictionary candidates; not validated NER",
             "clinical_efficacy": None, "species": None, "semantic_entailment_validated": False,
             "expert_review_status": "pending"} for row in relations]


def prepare(articles, source_manifest, ctd_dir, out, expected_documents=1198):
    out = Path(out)
    if (out / "preparation_lock.json").exists():
        raise FileExistsError("Literature preparation already frozen; choose a new versioned directory")
    manifest = json.loads(Path(source_manifest).read_text(encoding="utf-8"))
    if set(manifest["searches"]) != TOPICS or manifest["max_publication_date"] != "2026-10-01":
        raise ValueError("Expected exactly the six approved topics with October1 cutoff")
    if manifest.get("unique_documents") != expected_documents:
        raise ValueError("Frozen corpus count differs from specified expected_documents")
    if sha256(articles) != manifest["articles_sha256"]:
        raise ValueError("Articles do not match their frozen PubMed source manifest")
    snapshot = ctd_dictionary(ctd_dir)
    counts = Counter(record["type"] for record in snapshot["records"])
    if counts != {"CHEMICAL": 893, "GENE": 217}:
        raise ValueError(f"Expected scoped CTD 893chemical/217gene IDs, received {dict(counts)}")
    dictionary = MentionDictionary(snapshot["records"])
    documents = documents_from_articles(articles, dictionary)
    if len(documents) != expected_documents:
        raise ValueError("Actual frozen corpus record count mismatch")
    requested = {pmid for search in manifest["searches"].values() for pmid in search["selected_pmids"]}
    if {document["pmid"] for document in documents} != requested - set(manifest["missing_pmids"]):
        raise ValueError("Corpus PMID set does not match the frozen query lists/completeness ledger")
    out.mkdir(parents=True, exist_ok=True)
    atomic_json(out / "dictionary.json", snapshot)
    atomic_json(out / "dictionary_exclusions.json", {"excluded_aliases": dictionary.excluded, "shared_aliases": dictionary.ambiguous})
    with (out / "documents.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for document in documents:
            stream.write(json.dumps(document, ensure_ascii=False) + "\n")
    rules = [row for document in documents for row in annotate_relations(rule_extract(document), document)]
    atomic_json(out / "rule_relations.json", rules)
    write_relations(out / "rule_relations.tsv", rules)
    ledger = [{"pmid": document["pmid"], "eligible": document["qwen_eligible"],
               "abstention_reason": document["abstention_reason"], "entities": len(document["entities"]),
               "same_sentence_candidate_pairs": document["same_sentence_candidate_pairs"],
               "ambiguous_mentions": len(document["ambiguity_ledger"])} for document in documents]
    atomic_json(out / "eligibility_ledger.json", ledger)
    summary = {"documents": len(documents), "qwen_eligible_documents": sum(document["qwen_eligible"] for document in documents),
               "explicit_empty_abstention_documents": sum(not document["qwen_eligible"] for document in documents),
               "rule_relations": len(rules), "dictionary_ids": dict(counts),
               "mention_counts": dict(Counter(entity["type"] for document in documents for entity in document["entities"])),
               "ambiguity_mention_count": sum(len(document["ambiguity_ledger"]) for document in documents),
               "documents_with_same_sentence_candidate_pairs": sum(document["same_sentence_candidate_pairs"] > 0 for document in documents),
               "gold_available": False, "food_domain_precision": None, "ner_accuracy": None,
               "clinical_relations_extracted": False, "binding_affinity_labels": 0,
               "scope": "Scoped lexical chemical-gene candidates in six-query PubMed corpus; rule output requires expert review; no clinical efficacy or validated NER claim"}
    atomic_json(out / "preparation_summary.json", summary)
    lock = {"prepared_utc": now(), "articles_sha256": sha256(articles), "pubmed_source_manifest_sha256": sha256(source_manifest),
            "module_sha256": sha256(Path(__file__)), "expected_documents": expected_documents, "policy": summary["scope"],
            "gene_alias_policy": "Canonical CTD symbol (case-sensitive) and gene name only; unreviewed source synonyms retained/excluded; clinical terms blocked",
            "files": {name: sha256(out / name) for name in ("documents.jsonl", "dictionary.json", "dictionary_exclusions.json", "rule_relations.json", "rule_relations.tsv", "eligibility_ledger.json", "preparation_summary.json")},
            "source_rights": snapshot["rights"], "pubmed_abstract_rights": manifest["license_note"]}
    atomic_json(out / "preparation_lock.json", lock)
    return summary


def load_prepared(directory):
    directory = Path(directory)
    lock = json.loads((directory / "preparation_lock.json").read_text(encoding="utf-8"))
    for name, expected in lock["files"].items():
        if Path(name).name != name or sha256(directory / name) != expected:
            raise ValueError("Prepared literature artifact changed: " + name)
    # A literal Unicode U+2028/U+2029 in an original abstract is legal JSON
    # string content; str.splitlines would incorrectly split that source text.
    with (directory / "documents.jsonl").open(encoding="utf-8") as stream:
        documents = [json.loads(line) for line in stream]
    if len(documents) != lock["expected_documents"] or len({row["pmid"] for row in documents}) != len(documents):
        raise ValueError("Frozen prepared cohort count/PMIDs changed")
    for document in documents:
        if hashlib.sha256(document["text"].encode()).hexdigest() != document["text_sha256"]:
            raise ValueError("Original document text changed")
        for entity in document["entities"]:
            if document["text"][entity["start"]:entity["end"]] != entity["text"]:
                raise ValueError("Original candidate offsets changed")
    return documents, lock


def analyze_literature(documents, responses, out):
    out = Path(out)
    statuses, rejection_reasons = Counter(), Counter()
    reports, accepted = [], []
    for document in documents:
        row = responses.get(document["pmid"])
        if row is None:
            statuses["missing_response"] += 1
            reports.append({"pmid": document["pmid"], "accepted": [], "rejected": [], "error": "missing_response"})
            continue
        if not isinstance(row, dict) or not isinstance(row.get('status'), str):
            raise ValueError('Literature response requires an explicit status record')
        if not document["qwen_eligible"]:
            if row.get("status") != "no_comention_abstention" or row.get("response") != '{"relations": []}':
                raise ValueError("Noneligible document must remain an explicit empty abstention")
            statuses["no_comention_abstention"] += 1
            reports.append({"pmid": document["pmid"], "accepted": [], "rejected": [], "status": "no_comention_abstention"})
            continue
        statuses[row["status"]] += 1
        if row['status'] != 'generated':
            reports.append({'pmid': document['pmid'], 'accepted': [], 'rejected': [], 'error': row['status']})
            continue
        try:
            parsed = parse_qwen(row["response"], document)
            parsed["accepted"] = annotate_relations(parsed["accepted"], document)
            accepted.extend(parsed["accepted"])
            for rejection in parsed["rejected"]:
                rejection_reasons.update(rejection["reasons"])
        except (ValueError, TypeError, KeyError) as error:
            rejection_reasons["invalid_json_or_schema"] += 1
            parsed = {"pmid": document["pmid"], "accepted": [], "rejected": [], "error": str(error)}
        reports.append(parsed)
    atomic_json(out / "qwen_relations.json", accepted)
    write_relations(out / "qwen_relations.tsv", accepted)
    atomic_json(out / "validation_reports.json", reports)
    summary = {"documents": len(documents), "eligible_documents": sum(document["qwen_eligible"] for document in documents),
               "response_status_counts": dict(statuses), "accepted_relation_rows": len(accepted),
               "rejected_relation_rows": sum(len(row["rejected"]) for row in reports),
               "rejection_reason_counts": dict(rejection_reasons), "gold_available": False,
               "food_precision": None, "ner_accuracy": None, "semantic_entailment_validated": False,
               "clinical_efficacy_relations": 0, "experimental_binding_labels": 0,
               "interpretation": "Span-validated extracted candidates; dictionary mention and scientific entailment require expert review; public PubMed may be in pretraining"}
    atomic_json(out / "inference_summary.json", summary)
    return summary
