import pytest
from nutriomics_dti.protein_split_v2 import qualifying_edges,edge_components


def test_reported_nonrepresentative_edge_merges_entire_family(tmp_path):
    path = tmp_path/"alignments.tsv"
    path.write_text("A\tB\t0.61\t0.9\t0.95\nB\tC\t0.42\t0.81\t0.85\nC\tD\t0.99\t0.3\t0.9\n")
    edges = qualifying_edges(path)
    clusters = edge_components(["A","B","C","D"],edges)
    assert clusters["A"]==clusters["B"]==clusters["C"]
    assert clusters["D"]!=clusters["A"]


def test_malformed_fraction_or_foreign_target_is_not_silently_clustered(tmp_path):
    path = tmp_path/"bad.tsv"
    path.write_text("A\tB\t61\t0.9\t0.95\n")
    with pytest.raises(ValueError,match="fractions"):
        qualifying_edges(path)
    with pytest.raises(ValueError,match="outside"):
        edge_components(["A"],[("A","unknown")])
