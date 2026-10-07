import hashlib
import json
import pytest
from nutriomics_evidence.literature import MentionDictionary, documents_from_articles, analyze_literature, annotate_relations, load_prepared
from nutriomics_evidence.extraction import rule_extract
from nutriomics_evidence.sources import sha256, atomic_json


def lexicon():
    return MentionDictionary([
        {"id": "MESH:C1", "type": "CHEMICAL", "aliases": [{"text": "caffeine", "case_sensitive": False, "basis": "source chemical name"}]},
        {"id": "MESH:C2", "type": "CHEMICAL", "aliases": [{"text": "caffeine citrate", "case_sensitive": False, "basis": "source chemical name"}]},
        {"id": "NCBIGene:1", "type": "GENE", "aliases": [{"text": "ACE", "case_sensitive": True, "basis": "source gene symbol"}, {"text": "angiotensin converting enzyme", "case_sensitive": False, "basis": "source gene name"}]},
        {"id": "NCBIGene:2", "type": "GENE", "aliases": [{"text": "DBP", "case_sensitive": True, "basis": "ambiguous source symbol"}]},
    ])


def test_longest_offsets_case_sensitive_symbols_and_clinical_outcome_exclusion():
    text = 'Caffeine citrate inhibits ACE. Caffeine does not inhibit ace. DBP and SBP are blood pressure outcomes.'
    dictionary = lexicon()
    mentions, ledger = dictionary.annotate(text)
    assert [(entity["text"], entity["normalized_ids"]) for entity in mentions] == [
        ("Caffeine citrate", ["MESH:C2"]), ("ACE", ["NCBIGene:1"]), ("Caffeine", ["MESH:C1"])]
    assert all(text[entity["start"]:entity["end"]] == entity["text"] for entity in mentions)
    assert ledger == []
    assert dictionary.excluded[0]["reason"].startswith("clinical_outcome")


def test_shared_alias_stays_in_ambiguity_ledger_without_promoting_an_identity():
    dictionary = MentionDictionary([
        {"id": identifier, "type": "CHEMICAL", "aliases": [{"text": "vitamin D", "case_sensitive": False, "basis": "source synonym"}]}
        for identifier in ["MESH:D1", "MESH:D2"]])
    entities, ledger = dictionary.annotate("Vitamin D inhibits ACE.")
    assert entities == []
    assert ledger[0]["candidate_ids"] == ["MESH:D1", "MESH:D2"]


def test_native_gene_synonyms_the_diabetes_epa_do_not_become_gene_mentions():
    dictionary = MentionDictionary([
        {"id": "NCBIGene:7054", "type": "GENE", "aliases": [
            {"text": "TH", "case_sensitive": True, "basis": "CTD canonical gene symbol"},
            {"text": "The", "case_sensitive": True, "basis": "CTD gene synonym"}]},
        {"id": "NCBIGene:3953", "type": "GENE", "aliases": [{"text": "diabetes", "case_sensitive": False, "basis": "CTD gene synonym"}]},
        {"id": "NCBIGene:7076", "type": "GENE", "aliases": [{"text": "EPA", "case_sensitive": True, "basis": "CTD gene synonym"}]},
    ])
    entities, _ = dictionary.annotate("The EPA trial in diabetes lowered BP.")
    assert entities == []
    assert len([row for row in dictionary.excluded if row["reason"].startswith("unreviewed")]) == 3


def test_literal_title_space_abstract_offsets_rule_scope_and_explicit_abstention(tmp_path):
    articles = tmp_path / "articles.jsonl"
    rows = [{"pmid": "123", "title": "Mechanism.", "abstract": '"Caffeine" inhibits ACE.', "topics": ["caffeine"], "source_url": "https://pubmed.ncbi.nlm.nih.gov/123/"},
            {"pmid": "124", "title": "DASH and blood pressure.", "abstract": "DASH lowers SBP and DBP.", "topics": ["DASH"], "source_url": "https://pubmed.ncbi.nlm.nih.gov/124/"}]
    articles.write_text("\n".join(json.dumps(row) for row in rows), encoding="utf8")
    documents = documents_from_articles(articles, lexicon())
    first, second = documents
    assert first["text"] == 'Mechanism. "Caffeine" inhibits ACE.'
    assert first["qwen_eligible"] and not second["qwen_eligible"]
    relations = annotate_relations(rule_extract(first), first)
    assert len(relations) == 1 and relations[0]["subject_candidate_ids"] == ["MESH:C1"]
    assert relations[0]["clinical_efficacy"] is None and relations[0]["binding_affinity"] is None
    summary = analyze_literature(documents, {
        "123": {"status": "generated", "response": '{"relations": []}'},
        "124": {"status": "no_comention_abstention", "response": '{"relations": []}'},
    }, tmp_path)
    assert summary["documents"] == 2 and summary["response_status_counts"]["no_comention_abstention"] == 1
    assert summary["gold_available"] is False and summary["food_precision"] is None
    with pytest.raises(ValueError, match="explicit empty abstention"):
        analyze_literature(documents, {"124": {"status": "generated", "response": '{"relations": []}'}}, tmp_path)


def test_literal_unicode_paragraph_separator_does_not_split_jsonl_record(tmp_path):
    text = "Caffeine\u2028inhibits ACE."
    entities, _ = lexicon().annotate(text)
    document = {"pmid": "123", "text": text, "text_sha256": hashlib.sha256(text.encode()).hexdigest(), "entities": entities}
    path = tmp_path / "documents.jsonl"
    path.write_text(json.dumps(document, ensure_ascii=False) + "\n", encoding="utf8")
    atomic_json(tmp_path / "preparation_lock.json", {"expected_documents": 1, "files": {"documents.jsonl": sha256(path)}})
    documents, _ = load_prepared(tmp_path)
    assert len(documents) == 1 and documents[0]["text"] == text
    assert documents[0]["entities"][1]["text"] == "ACE"
