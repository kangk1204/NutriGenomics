import pandas as pd
import pytest
from nutriomics_dti.splits import assert_no_overlap,protein_clusters,kmers,kmer_jaccard


def frame(pair,source,compound="c",target="t",cluster="cl"):
    return pd.DataFrame([{"pair_id":pair,"source_group":source,"compound_id":compound,"target_id":target,"protein_cluster":cluster}])


def test_source_and_pair_overlap_are_rejected():
    train,val,test = frame("a","paper1"),frame("b","paper2"),frame("c","paper1")
    with pytest.raises(ValueError,match="source_group"):
        assert_no_overlap(train,val,test,"random")
    test = frame("a","paper3")
    with pytest.raises(ValueError,match="pair_id"):
        assert_no_overlap(train,val,test,"random")


def test_similarity_groups_are_transitive_and_cannot_cross_folds():
    seq = "ACDEFGHIKLMNPQRSTVWY"*3
    proteins = pd.DataFrame({"target_id":["a","b","c"],"protein_sequence":[seq,seq[:-1]+"A","YYYYYYYYYYYYYYYYYYYY"]})
    groups,edges = protein_clusters(proteins,0.5)
    assert groups["a"]==groups["b"]
    assert groups["a"]!=groups["c"]
    assert edges>=1
    assert kmer_jaccard(kmers(seq),kmers(seq))==1
    with pytest.raises(ValueError,match="protein_cluster"):
        assert_no_overlap(frame("a","x",target="a",cluster="same"),frame("b","y",target="b",cluster="other"),frame("c","z",target="c",cluster="same"),"protein_similarity")
