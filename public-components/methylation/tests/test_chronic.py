import numpy as np
import pytest
from nutriomics_methylation.chronic import TrainingProbeSelector,nested,tune


def test_chronic_training_filter_ignores_holdout_distribution():
    X=np.array([[.1,.2,np.nan,1],[.2,.4,.2,1],[.4,.3,.3,1],[.8,.7,.4,1],[.7,.6,.5,1],[.9,.9,.6,1]])
    y=np.array([0,0,0,1,1,1])
    fitted=TrainingProbeSelector(k=10,max_missing=.1).fit(X,y)
    assert set(fitted.selected_indices_)=={0,1}
    med=fitted.medians_.copy();selected=fitted.selected_indices_.copy()
    transformed=fitted.transform(np.array([[np.nan,999,-999,0]]))
    np.testing.assert_array_equal(med,fitted.medians_)
    np.testing.assert_array_equal(selected,fitted.selected_indices_)
    assert transformed[0,np.where(selected==0)[0][0]]==pytest.approx(np.median(X[:,0]),rel=1e-6)
    assert fitted.n_fit_samples_==6
    with pytest.raises(ValueError,match='width'):fitted.transform(np.zeros((2,3)))


def test_chronic_nested_membership_and_inner_learned_denominators():
    X=np.random.default_rng(8).normal(size=(20,8));y=np.tile([0,1],10)
    ids=np.array([f's{i}' for i in range(20)])
    p,fold,audits=nested(X,y,ids,17,ranking='variance')
    assert np.isfinite(p).all() and set(fold)==set(range(5))
    test_ids=[]
    for a in audits:
        assert not set(a['train_sample_ids'])&set(a['test_sample_ids'])
        assert a['preprocessing_fit_n']==len(a['train_sample_ids'])
        test_ids+=a['test_sample_ids']
        for inner in a['inner_splits']:
            assert not set(inner['train_indices'])&set(inner['validation_indices'])
            assert inner['preprocessing_fit_n']==len(inner['train_indices'])
    assert sorted(test_ids)==sorted(ids)


def test_chronic_refuses_invalid_outcomes_and_empty_training_features():
    with pytest.raises(ValueError,match='binary'):TrainingProbeSelector().fit(np.ones((3,2)),[0,0,0])
    with pytest.raises(ValueError,match='No usable'):TrainingProbeSelector().fit(np.ones((4,2)),[0,0,1,1])
    with pytest.raises(ValueError,match='Insufficient'):tune(np.arange(12).reshape(4,3),np.array([0,0,0,1]),1)


def test_duplicate_sample_identifiers_cannot_cross_chronic_folds():
    X=np.random.default_rng(10).normal(size=(20,8));y=np.tile([0,1],10)
    with pytest.raises(ValueError,match='overlap'):nested(X,y,np.array(['same']*20),17)
