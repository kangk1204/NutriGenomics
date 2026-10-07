import numpy as np
import pytest
from nutriomics_methylation.extension_grid import splits,tune_all,calibrated_fit,nested_all

def dataset():
    rng=np.random.default_rng(37);y=np.tile([0,1],18);X=rng.normal(size=(36,24)).astype(np.float32)
    X[:,0]+=y*2;X[0,2]=np.nan;return X,y

def test_grouped_folds_never_share_family():
    X,y=dataset();groups=np.repeat(np.arange(18),2)
    for train,test in splits(y,11,5,groups):assert not set(groups[train])&set(groups[test])

def test_calibration_preprocessing_is_inside_score_training():
    X,y=dataset();choices=tune_all(X,y,11,ks=(2,5),Cs=(.1,1.))
    for family,choice in choices.items():
        fitted=calibrated_fit(X,y,choice,11)
        for fold in fitted.audit['calibration_splits']:
            assert set(fold['train_indices']).isdisjoint(fold['calibration_score_indices'])
            assert fold['preprocessing_fit_n']==len(fold['train_indices'])
        assert np.allclose(fitted.predict_proba(X).sum(axis=1),1)

def test_outer_predictions_and_masks_complete():
    X,y=dataset();samples=np.array([f'P{i}' for i in range(len(y))])
    result=nested_all(X,y,samples,11,ks=(2,),Cs=(.1,))
    for item in result.values():
        assert np.isfinite(item['probability']).all()
        for audit in item['audits']:
            assert set(audit['outer_train_samples']).isdisjoint(audit['outer_test_samples'])
            assert audit['calibration']['base_fit_n']==len(audit['outer_train_samples'])

def test_duplicate_patient_ids_fail():
    X,y=dataset()
    with pytest.raises(ValueError,match='Repeated patient'):nested_all(X,y,np.repeat('same',len(y)),11,ks=(2,),Cs=(.1,))

def test_confounding_center_partition_is_feasible_without_subject_fallback():
    counts={1:(0,26),3:(18,0),4:(3,5),5:(0,4),7:(21,2),8:(0,18)}
    y=np.concatenate([np.repeat([0,1],c) for c in counts.values()]);groups=np.concatenate([np.repeat(g,sum(c)) for g,c in counts.items()])
    outer=splits(y,20261003,5,groups)
    assert len(outer)==3
    for train,test in outer:
        assert set(y[train])==set(y[test])=={0,1}
        assert not set(groups[train])&set(groups[test])
        assert len(splits(y[train],20261003,3,groups[train]))==2
