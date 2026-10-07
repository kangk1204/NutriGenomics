import pytest
from nutriomics_evidence.token_evaluation import spans,span_metrics,encode_word_chunks,reconstruct_word_predictions


def test_exact_boundaries_and_type_changes():
    assert spans(['I-GENE','I-GENE','B-GENE','O'])=={(0,2,'GENE'),(2,3,'GENE')}
    metrics=span_metrics([['B-GENE','I-GENE','O']],[['B-GENE','O','O']])
    assert metrics['precision']==0 and metrics['false_positive_spans']==1 and metrics['false_negative_spans']==1
    with pytest.raises(ValueError):span_metrics([['O']],[[]])


class Encoding(dict):
    def __init__(self,word_ids):
        super().__init__(input_ids=list(range(len(word_ids))),attention_mask=[1]*len(word_ids));self.ids=word_ids
    def word_ids(self):return self.ids


class FakeTokenizer:
    is_fast=True
    def num_special_tokens_to_add(self,pair=False):return 2
    def __call__(self,words,add_special_tokens=True,**kwargs):
        ids=[]
        for i,word in enumerate(words):ids.extend([i]*(2 if word=='splitgene' else 1))
        return Encoding(([None]+ids+[None]) if add_special_tokens else ids)


def test_long_sentence_cross_chunk_gene_and_first_subword_reconstruction():
    rows=[{'words':['plain','splitgene','next','tail'],'labels':['O','B-GENE','I-GENE','O']}]
    pieces=encode_word_chunks(FakeTokenizer(),rows,['O','B-GENE','I-GENE'],True,max_tokens=5)
    assert len(pieces)==2 and pieces[0]['features']['labels']==[-100,0,1,-100,-100]
    tags=[[0,0,1,0,0],[0,2,0,0]]
    reconstructed=reconstruct_word_predictions(rows,pieces,tags,['O','B-GENE','I-GENE'])
    assert reconstructed==[rows[0]['labels']]
    assert spans(reconstructed[0])=={(1,3,'GENE')}


def test_truncated_or_missing_word_prediction_and_invalid_IOB_are_rejected():
    rows=[{'words':['x'],'labels':['O']}]
    pieces=encode_word_chunks(FakeTokenizer(),rows,['O','B-GENE','I-GENE'],False,max_tokens=5)
    with pytest.raises(ValueError):reconstruct_word_predictions(rows,pieces,[[0]],['O','B-GENE','I-GENE'])
    with pytest.raises(ValueError):reconstruct_word_predictions(rows,[],[],['O'])
    with pytest.raises(ValueError):spans(['BAD-GENE'])
