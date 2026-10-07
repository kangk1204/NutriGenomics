"""Exploratory chronic-disease screening with fully nested learned preprocessing."""
from __future__ import annotations
import warnings
import numpy as np
from sklearn.base import BaseEstimator,TransformerMixin
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class TrainingProbeSelector(TransformerMixin,BaseEstimator):
    def __init__(self,k=200,max_missing=0.1,ranking='supervised'):
        self.k=k;self.max_missing=max_missing;self.ranking=ranking

    def fit(self,X,y):
        X=np.asarray(X,dtype=np.float32);y=np.asarray(y)
        if X.ndim!=2 or len(X)!=len(y) or set(y)!={0,1}:raise ValueError('Aligned binary training data required')
        if self.ranking not in ['supervised','variance']:raise ValueError('Unknown ranking')
        self.n_fit_samples_=len(X);self.n_features_in_=X.shape[1]
        keep=np.flatnonzero(np.mean(~np.isfinite(X),axis=0)<=self.max_missing)
        clean=np.where(np.isfinite(X[:,keep]),X[:,keep],np.nan)
        med=np.nanmedian(clean,axis=0)
        clean=np.where(np.isfinite(clean),clean,med)
        var=np.var(clean,axis=0,dtype=np.float64);nonconstant=var>1e-12
        keep=keep[nonconstant];med=med[nonconstant];clean=clean[:,nonconstant];var=var[nonconstant]
        if not len(keep):raise ValueError('No usable training probes')
        if self.ranking=='supervised':
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                scores=f_classif(clean,y)[0]
            scores=np.nan_to_num(scores,nan=-np.inf,posinf=np.finfo(float).max)
        else:scores=var
        order=np.argsort(-scores,kind='stable')[:min(self.k,len(keep))]
        self.selected_indices_=keep[order];self.medians_=med[order]
        return self

    def transform(self,X):
        X=np.asarray(X,dtype=np.float32)
        if X.shape[1]!=self.n_features_in_:raise ValueError('Probe width changed')
        values=X[:,self.selected_indices_]
        return np.where(np.isfinite(values),values,self.medians_)


def model(C,seed):
    return LogisticRegression(C=C,penalty='l2',solver='liblinear',class_weight='balanced',max_iter=3000,tol=1e-4,random_state=seed)


def tune(X,y,seed,ranking='supervised',ks=(10,50,200),Cs=(0.1,1.0)):
    nfolds=min(3,int(np.bincount(y).min()))
    if nfolds<2:raise ValueError('Insufficient minority participants for inner CV')
    grid=[(k,C) for k in ks for C in Cs];scores=[[] for _ in grid];splits=[]
    for tr,va in StratifiedKFold(nfolds,shuffle=True,random_state=seed).split(X,y):
        selector=TrainingProbeSelector(max(ks),ranking=ranking).fit(X[tr],y[tr])
        a=selector.transform(X[tr]);b=selector.transform(X[va])
        for k in ks:
            width=min(k,a.shape[1]);scaler=StandardScaler().fit(a[:,:width])
            train=scaler.transform(a[:,:width]);valid=scaler.transform(b[:,:width])
            for C in Cs:
                fitted=model(C,seed).fit(train,y[tr])
                scores[grid.index((k,C))].append(float(roc_auc_score(y[va],fitted.predict_proba(valid)[:,1])))
        splits.append({'train_indices':tr.tolist(),'validation_indices':va.tolist(),'preprocessing_fit_n':selector.n_fit_samples_})
    best=int(np.argmax([np.mean(s) for s in scores]));k,C=grid[best]
    pipe=Pipeline([('probes',TrainingProbeSelector(k,ranking=ranking)),('scale',StandardScaler()),('model',model(C,seed))]).fit(X,y)
    return pipe,{'k':k,'C':C,'ranking':ranking,'inner_folds':nfolds,'inner_splits':splits,'candidate_scores':[{'k':k0,'C':C0,'auroc':s} for (k0,C0),s in zip(grid,scores)]}


def nested(X,y,samples,seed,ranking='supervised'):
    count=min(5,int(np.bincount(y).min()));pred=np.full(len(y),np.nan);fold=np.full(len(y),-1);audits=[]
    for i,(tr,te) in enumerate(StratifiedKFold(count,shuffle=True,random_state=seed).split(X,y)):
        fitted,audit=tune(X[tr],y[tr],seed+i,ranking)
        if set(samples[tr])&set(samples[te]):raise ValueError('Participant overlap')
        if np.isfinite(pred[te]).any():raise ValueError('Repeated OOF assignment')
        pred[te]=fitted.predict_proba(X[te])[:,1];fold[te]=i
        audit.update({'outer_fold':i,'train_sample_ids':samples[tr].tolist(),'test_sample_ids':samples[te].tolist(),
                      'selected_indices':fitted['probes'].selected_indices_.tolist(),'preprocessing_fit_n':fitted['probes'].n_fit_samples_,
                      'all_learned_preprocessing_training_only':True})
        audits.append(audit)
    if not np.isfinite(pred).all():raise ValueError('Missing OOF predictions')
    return pred,fold,audits
