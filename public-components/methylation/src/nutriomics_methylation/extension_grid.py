"""Source-only nested classification with fold-local probe QC and calibration.

This module deliberately has no external outcome argument in tuning. All
reported probabilities are calibrated with out-of-fold scores within the
corresponding development population. No global variance/ANOVA mask is learned.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import product
import warnings
import numpy as np
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from .chronic import TrainingProbeSelector

FAMILIES=('ridge','elasticnet','linear_svm')
FEATURE_COUNTS=(10,50,200,1000)
REGULARIZATION=(.01,.1,1.,10.)

def estimator(family,C,seed):
    if family=='ridge':return LogisticRegression(C=C,penalty='l2',solver='liblinear',class_weight='balanced',max_iter=5000,random_state=seed)
    if family=='elasticnet':return LogisticRegression(C=C,penalty='elasticnet',l1_ratio=.5,solver='saga',class_weight='balanced',max_iter=5000,tol=1e-3,random_state=seed)
    if family=='linear_svm':return LinearSVC(C=C,class_weight='balanced',dual='auto',max_iter=10000,random_state=seed)
    raise ValueError('Unknown classifier family')

def splits(y,seed,n=3,groups=None):
    y=np.asarray(y,int)
    n=min(n,int(np.bincount(y).min()))
    if groups is not None:
        groups=np.asarray(groups)
        n=min(n,*[len(np.unique(groups[y==c])) for c in (0,1)])
    if n<2:raise ValueError('Insufficient independent groups per class')
    if groups is not None and len(np.unique(groups))<=10:
        # Small, strongly confounded recruitment-center designs need an exact
        # label-count partition. Exhaustively choose the most balanced feasible
        # assignment, with lexicographic tie resolution. No methylation or
        # prediction values enter this partition choice.
        identifiers=np.unique(groups);counts=np.array([[sum((groups==g)&(y==c)) for c in (0,1)] for g in identifiers]);best=None
        for rest in product(range(n),repeat=len(identifiers)-1):
            assignment=np.array((0,*rest));table=np.array([counts[assignment==f].sum(axis=0) for f in range(n)])
            if (table==0).any():continue
            objective=float(np.sum((table/counts.sum(axis=0)-1/n)**2))
            if best is None or objective<best[0]-1e-12:best=(objective,assignment)
        if best is None:raise ValueError('No feasible center-grouped partition')
        result=[]
        for f in range(n):
            test=np.flatnonzero(np.isin(groups,identifiers[best[1]==f]));train=np.flatnonzero(~np.isin(groups,identifiers[best[1]==f]));result.append((train,test))
        return result
    # Recruitment centers can be strongly aliased with subtype. Choose the first
    # feasible partition from a frozen sequence using labels/group counts only,
    # never features or predictive results. Save actual memberships in each audit.
    for offset in range(1 if groups is None else 100):
        cv=StratifiedKFold(n,shuffle=True,random_state=seed) if groups is None else StratifiedGroupKFold(n,shuffle=True,random_state=seed+offset)
        result=list(cv.split(np.zeros(len(y)),y,groups))
        if any(set(y[tr])!={0,1} or set(y[te])!={0,1} for tr,te in result):continue
        for train,test in result:
            if set(train)&set(test):raise ValueError('Index overlap')
            if groups is not None and set(groups[train])&set(groups[test]):raise ValueError('Group overlap')
        return result
    raise ValueError('No feasible grouped partition with both classes in every fold; no participant-fold fallback is allowed')

def pipe(family,k,C,seed):
    return Pipeline([('probes',TrainingProbeSelector(k,ranking='supervised')),('scale',StandardScaler()),('model',estimator(family,C,seed))])

def fit_checked(model,X,y):
    with warnings.catch_warnings(record=True) as observed:
        warnings.simplefilter('always',ConvergenceWarning)
        model.fit(X,y)
    return sum(issubclass(w.category,ConvergenceWarning) for w in observed)

def tune_all(X,y,seed,groups=None,ks=FEATURE_COUNTS,Cs=REGULARIZATION):
    grid=[(k,C) for k in ks for C in Cs]
    scores={f:[[] for _ in grid] for f in FAMILIES};counts={f:0 for f in FAMILIES};audit=[]
    for train,valid in splits(y,seed,3,groups):
        selector=TrainingProbeSelector(max(ks),ranking='supervised').fit(X[train],y[train])
        a,b=selector.transform(X[train]),selector.transform(X[valid])
        for k in ks:
            width=min(k,a.shape[1]);scaler=StandardScaler().fit(a[:,:width]);aa=scaler.transform(a[:,:width]);bb=scaler.transform(b[:,:width])
            for C in Cs:
                for family in FAMILIES:
                    fitted=estimator(family,C,seed);counts[family]+=fit_checked(fitted,aa,y[train])
                    score=float(roc_auc_score(y[valid],fitted.decision_function(bb)))
                    scores[family][grid.index((k,C))].append(score)
        audit.append({'train_indices':train.tolist(),'validation_indices':valid.tolist(),'preprocessing_fit_n':selector.n_fit_samples_,
                      'training_usable_probes':len(selector.selected_indices_),'group_disjoint':True})
    selections={}
    for family in FAMILIES:
        means=[np.mean(v) for v in scores[family]];best=int(np.argmax(means));k,C=grid[best]
        selections[family]={'family':family,'k':int(k),'C':float(C),'elasticnet_l1_ratio':.5 if family=='elasticnet' else None,
          'inner_mean_auroc':float(means[best]),'convergence_warning_count':counts[family],
          'candidate_scores':[{'k':int(k0),'C':float(c0),'fold_auroc':v,'mean_auroc':float(np.mean(v))} for (k0,c0),v in zip(grid,scores[family])],
          'inner_splits':audit,'tie_break':'Smallest k then smallest C; frozen grid order'}
    return selections

@dataclass
class CalibratedSourceModel:
    base: Pipeline
    calibrator: LogisticRegression
    audit: dict

    def predict_proba(self,X):
        z=np.asarray(self.base.decision_function(X)).reshape(-1,1)
        return self.calibrator.predict_proba(z)

def calibrated_fit(X,y,choice,seed,groups=None):
    score=np.full(len(y),np.nan);audits=[];warnings_total=0
    family,k,C=choice['family'],choice['k'],choice['C']
    for train,valid in splits(y,seed+107,3,groups):
        fitted=pipe(family,k,C,seed);warnings_total+=fit_checked(fitted,X[train],y[train])
        score[valid]=fitted.decision_function(X[valid])
        audits.append({'train_indices':train.tolist(),'calibration_score_indices':valid.tolist(),
                      'preprocessing_fit_n':fitted['probes'].n_fit_samples_,
                      'selected_indices':fitted['probes'].selected_indices_.tolist(),'group_disjoint':True})
    if not np.isfinite(score).all():raise ValueError('Incomplete calibration OOF')
    # Unweighted logistic mapping targets the source empirical prevalence.
    calibration=LogisticRegression(C=1e6,solver='lbfgs',max_iter=3000).fit(score.reshape(-1,1),y)
    fitted=pipe(family,k,C,seed);warnings_total+=fit_checked(fitted,X,y)
    return CalibratedSourceModel(fitted,calibration,{'calibration_method':f'{len(audits)}-fold OOF logistic score mapping; no external fitting',
      'calibration_fold_count':len(audits),
      'calibration_splits':audits,'calibrator_fit_n':len(y),'base_fit_n':len(y),'convergence_warning_count':warnings_total,
      'calibrator_coefficient':float(calibration.coef_[0,0]),'calibrator_intercept':float(calibration.intercept_[0])})

def nested_all(X,y,samples,seed,groups=None,ks=FEATURE_COUNTS,Cs=REGULARIZATION):
    X=np.asarray(X,np.float32);y=np.asarray(y,int);samples=np.asarray(samples)
    if len(np.unique(samples))!=len(samples):raise ValueError('Repeated patient IDs must be represented using an explicit group design')
    result={f:{'probability':np.full(len(y),np.nan),'fold':np.full(len(y),-1),'audits':[]} for f in FAMILIES}
    for index,(train,test) in enumerate(splits(y,seed,5,groups)):
        if set(samples[train])&set(samples[test]):raise ValueError('Patient leakage')
        train_groups=None if groups is None else np.asarray(groups)[train]
        choices=tune_all(X[train],y[train],seed+index,train_groups,ks,Cs)
        for family in FAMILIES:
            fitted=calibrated_fit(X[train],y[train],choices[family],seed+index,train_groups)
            pred=fitted.predict_proba(X[test])[:,1]
            item=result[family];item['probability'][test]=pred;item['fold'][test]=index
            item['audits'].append({'outer_fold':index,'outer_train_samples':samples[train].tolist(),'outer_test_samples':samples[test].tolist(),
               'outer_train_indices':train.tolist(),'outer_test_indices':test.tolist(),
               'choice':choices[family],'calibration':fitted.audit,'selected_indices':fitted.base['probes'].selected_indices_.tolist(),
               'all_learned_preprocessing_and_calibration_training_only':True})
        print(f'completed outer fold {index+1} seed {seed}',flush=True)
    for item in result.values():
        if not np.isfinite(item['probability']).all():raise ValueError('Missing source OOF')
    return result
