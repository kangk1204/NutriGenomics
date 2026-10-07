"""Native literal TSV quoting and strict original-offset regression tests."""
from pathlib import Path
import pytest
from nutriomics_evidence.extraction import load_drugprot


def write_fixture(tmp_path: Path, entity_start: int = 9):
    abstracts = tmp_path / "abstracts.tsv"
    entities = tmp_path / "entities.tsv"
    # Literal quotes occur at both title and abstract starts in official v1.2.
    title, abstract = '"Title".', '"Ecstasy" inhibits AKT1.'
    abstracts.write_text("123\t" + title + "\t" + abstract + "\n", encoding="utf-8")
    entities.write_text(
        f"123\tT1\tCHEMICAL\t{entity_start}\t{entity_start+9}\t\"Ecstasy\"\n"
        "123\tT2\tGENE\t28\t32\tAKT1\n", encoding="utf-8"
    )
    return abstracts, entities


def test_native_literal_quotes_preserve_original_entity_offsets(tmp_path):
    documents = load_drugprot(*write_fixture(tmp_path))
    assert documents[0]["text"] == '"Title". "Ecstasy" inhibits AKT1.'
    assert documents[0]["entities"][0]["text"] == '"Ecstasy"'
    assert documents[0]["text"][9:18] == '"Ecstasy"'
    assert documents[0]["text"][28:32] == "AKT1"


def test_native_offset_mismatch_is_rejected_without_repair(tmp_path):
    with pytest.raises(ValueError, match="entity offset/type/text mismatch"):
        load_drugprot(*write_fixture(tmp_path, entity_start=10))
