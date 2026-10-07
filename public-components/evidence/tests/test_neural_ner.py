import pytest
from nutriomics_evidence.neural_ner import native_word_offsets,mentions_from_labels


def test_unicode_source_offsets_and_multiword_gene_span():
    text='α  tumor necrosis factor (TNF-alpha).'
    tokens=native_word_offsets(text)
    labels=['O']*len(tokens)
    for index,token in enumerate(tokens):
        if token['word']=='tumor':labels[index]='B-GENE'
        if token['word'] in ['necrosis','factor']:labels[index]='I-GENE'
    mentions=mentions_from_labels(text,tokens,labels)
    assert len(mentions)==1 and mentions[0]['text']=='tumor necrosis factor'
    assert text[mentions[0]['start']:mentions[0]['end']]==mentions[0]['text']


def test_offset_mismatch_and_empty_source():
    assert native_word_offsets(' \n')==[]
    with pytest.raises(ValueError):mentions_from_labels('ACE',[{'word':'ACE','start':1,'end':3}],['B-GENE'])


def test_astral_unicode_offsets_are_codepoints_not_utf16_or_bytes():
    text='🧬 ACE';tokens=native_word_offsets(text)
    mention=mentions_from_labels(text,tokens,['O','B-GENE'])[0]
    assert mention['start']==2 and mention['end']==5 and mention['text']=='ACE'
    assert mention['offset_unit']=='Unicode codepoints'
