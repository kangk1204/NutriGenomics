from pathlib import Path
import pytest
from nutriomics_evidence.neural_data import mark_context, pair_rows, read_cpr_gold


def test_markers_preserve_source_order_and_both_entities():
    text = 'ACE is inhibited by quercetin.'
    gene = {'start':0,'end':3,'text':'ACE'}
    chemical = {'start':20,'end':29,'text':'quercetin'}
    # Fix literal coordinates against source, not a shifted token string.
    chemical['start'] = text.index('quercetin')
    chemical['end'] = chemical['start']+len(chemical['text'])
    rendered = mark_context(text,chemical,gene)
    assert '[GENE] ACE [/GENE]' in rendered
    assert '[CHEM] quercetin [/CHEM]' in rendered
    assert 'CPR:' not in rendered
    bad = dict(gene,end=4)
    with pytest.raises(ValueError,match='offset'):
        mark_context(text,chemical,bad)


def test_native_multilabel_relations_kept(tmp_path):
    path = tmp_path/'gold.tsv'
    path.write_text('1\tCPR:3\tArg1:T1\tArg2:T2\n1\tCPR:4\tArg1:T1\tArg2:T2\n')
    assert read_cpr_gold(path)[('1','T1','T2')] == {'CPR:3','CPR:4'}
    path.write_text('1\tCPR:1\tArg1:T1\tArg2:T2\n')
    with pytest.raises(ValueError,match='native'):
        read_cpr_gold(path)
