import platform
import random
import time
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import ExtraTreesRegressor
from .features import features, applicability
from .io import load_json, sha256, utc_now, write_json
from .splits import load_locked,alignment_audit


def metrics(actual,predicted):
    actual,predicted = np.asarray(actual),np.asarray(predicted)
    if actual.ndim != 1 or predicted.ndim != 1 or actual.shape != predicted.shape or actual.size == 0:
        raise ValueError("Metrics require nonempty aligned one-dimensional labels and predictions")
    valid = np.isfinite(actual)&np.isfinite(predicted)
    if not valid.all():
        raise ValueError("Nonfinite prediction/label")
    correlation = float(spearmanr(actual,predicted).statistic) if len(actual)>1 and np.std(predicted)>0 and np.std(actual)>0 else None
    return {"n":len(actual),"rmse_pkd":float(np.sqrt(np.mean((actual-predicted)**2))),
            "mae_pkd":float(np.mean(np.abs(actual-predicted))),"spearman":correlation}


def train(split,out,model_name="baseline",seed=17,epochs=40,patience=6,device="auto",trees=150):
    (training,validation,_),lock = load_locked(split)
    out = Path(out)
    if (out/"training.json").exists():
        raise FileExistsError(f"Model already trained: {out}; preserve frozen run")
    out.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    history = []
    runtime = {"python":platform.python_version(),"platform":platform.platform(),"numpy":np.__version__}
    if model_name == "baseline":
        model = ExtraTreesRegressor(n_estimators=trees,max_depth=24,min_samples_leaf=2,max_features=0.7,
                                    random_state=seed,n_jobs=6)
        model.fit(features(training),training.pkd)
        predicted = model.predict(features(validation))
        joblib.dump(model,out/"baseline.joblib")
        import sklearn
        runtime["sklearn"] = sklearn.__version__
        config = {"trees":trees,"max_depth":24,"min_samples_leaf":2,"max_features":0.7,
                  "features":"Morgan radius2 chiral 1024bits +20AA composition+400dipeptides+log-length"}
    elif model_name == "cnn":
        import torch
        from .neural import AffinityCNN,PairDataset,predict_neural
        random.seed(seed)
        np.random.seed(seed)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.benchmark = False
        torch.use_deterministic_algorithms(True,warn_only=True)
        torch.set_num_threads(6)
        selected = "cuda" if device=="auto" and torch.cuda.is_available() else ("cpu" if device=="auto" else device)
        if selected=="cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA requested but unavailable")
        model = AffinityCNN().to(selected)
        mean = float(training.pkd.mean())
        scale = max(float(training.pkd.std()),1e-6)
        loader = torch.utils.data.DataLoader(PairDataset(training,mean,scale),batch_size=128,shuffle=True,
                                            generator=torch.Generator().manual_seed(seed),num_workers=0)
        optimizer = torch.optim.AdamW(model.parameters(),lr=1e-3,weight_decay=1e-4)
        best_loss = float("inf")
        stale = 0
        best_state = None
        for epoch in range(1,epochs+1):
            model.train()
            losses = []
            for drug,protein,label in loader:
                optimizer.zero_grad(set_to_none=True)
                prediction = model(drug.to(selected),protein.to(selected))
                loss = torch.nn.functional.mse_loss(prediction,label.to(selected))
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(),5)
                optimizer.step()
                losses.append(float(loss.detach().cpu()))
            predicted = predict_neural(model,validation,selected,mean,scale)
            score = metrics(validation.pkd,predicted)
            history.append({"epoch":epoch,"train_normalized_mse":float(np.mean(losses)),"validation":score})
            print(f"epoch={epoch} val_rmse={score['rmse_pkd']:.4f}",flush=True)
            if score["rmse_pkd"] < best_loss:
                best_loss = score["rmse_pkd"]
                best_state = {key:value.detach().cpu().clone() for key,value in model.state_dict().items()}
                stale = 0
            else:
                stale += 1
                if stale>=patience:
                    break
        model.load_state_dict(best_state)
        torch.save({"state_dict":best_state,"mean":mean,"scale":scale},out/"cnn.pt")
        predicted = predict_neural(model,validation,selected,mean,scale)
        runtime.update({"torch":torch.__version__,"cuda_build":torch.version.cuda,"device":selected,
                        "gpu":torch.cuda.get_device_name(0) if selected=="cuda" else None})
        config = {"fresh_initialization":True,"epochs_max":epochs,"patience":patience,"batch_size":128,
                  "drug_max_characters":256,"protein_max_residues":1024,"learning_rate":0.001,
                  "sequence_policy":"first 1024 residues; full sequence retained for identity and applicability checks",
                  "normalization_train_mean":mean,"normalization_train_std":scale}
    else:
        raise ValueError(model_name)
    training[["compound_id","target_id","smiles","protein_sequence","pair_id"]].drop_duplicates().to_csv(out/"training_entities.csv.gz",index=False)
    record = {"trained_utc":utc_now(),"elapsed_seconds":time.monotonic()-started,"model":model_name,"model_seed":seed,
              "split_mode":lock["mode"],"split_seed":lock["seed"],"split_lock_sha256":sha256(Path(split)/"lock.json"),
              "locked_test_sha256":lock["files"]["test"]["sha256"],"validation_metrics":metrics(validation.pkd,predicted),
              "config":config,"runtime":runtime,"history":history,"model_selection":"fixed architecture; CNN best validation epoch only; no test tuning"}
    write_json(out/"training.json",record)
    return record


def predict(model_dir,frame,device="auto"):
    model_dir = Path(model_dir)
    record = load_json(model_dir/"training.json")
    if record["model"]=="baseline":
        prediction = joblib.load(model_dir/"baseline.joblib").predict(features(frame))
    else:
        import torch
        from .neural import AffinityCNN,predict_neural
        selected = "cuda" if device=="auto" and torch.cuda.is_available() else ("cpu" if device=="auto" else device)
        checkpoint = torch.load(model_dir/"cnn.pt",map_location="cpu",weights_only=True)
        model = AffinityCNN().to(selected)
        model.load_state_dict(checkpoint["state_dict"])
        prediction = predict_neural(model,frame,selected,checkpoint["mean"],checkpoint["scale"])
    training = pd.read_csv(model_dir/"training_entities.csv.gz",keep_default_na=False)
    output = applicability(training,frame)
    output["predicted_pkd"] = prediction
    output["predicted_kd_nm"] = np.power(10.0,9-prediction)
    output["evidence_status"] = "model prediction; not an experimentally measured interaction"
    return output


def evaluate(split,model_dir,out,device="auto"):
    (training,_,test),lock = load_locked(split)
    record = load_json(Path(model_dir)/"training.json")
    if record["split_lock_sha256"]!=sha256(Path(split)/"lock.json"):
        raise ValueError("Model trained against a different locked split")
    output = predict(model_dir,test,device)
    output["observed_pkd"] = test.pkd
    output["source_group"] = test.source_group
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    output.to_csv(out/"test_predictions.csv.gz",index=False)
    result = {"evaluated_utc":utc_now(),"model":record["model"],"model_seed":record["model_seed"],
              "split_mode":lock["mode"],"split_seed":lock["seed"],"test_sha256":lock["files"]["test"]["sha256"],
              "test_metrics":metrics(test.pkd,output.predicted_pkd),"overlap_audit":lock["overlap_counts"],
              "applicability_fraction":float(output.applicability_flag.mean()),
              "max_train_morgan_tanimoto_quantiles":{str(q):float(output.max_train_morgan_tanimoto.quantile(q)) for q in (0,.25,.5,.75,1)},
              "max_train_protein_3mer_jaccard_quantiles":{str(q):float(output.max_train_protein_3mer_jaccard.quantile(q)) for q in (0,.25,.5,.75,1)},
              "evaluation_type":"internal provenance-purged locked holdout; not independent external validation"}
    ranks = []
    for _, group in output.groupby("target_id"):
        if group.compound_id.nunique()>=5 and group.observed_pkd.nunique()>1 and group.predicted_pkd.nunique()>1:
            ranks.append(float(spearmanr(group.observed_pkd,group.predicted_pkd).statistic))
    result["within_target_ranking"] = {"minimum_distinct_compounds":5,"eligible_targets":len(ranks),
                                        "macro_spearman_mean":float(np.mean(ranks)) if ranks else None,
                                        "macro_spearman_median":float(np.median(ranks)) if ranks else None}
    if lock["mode"]=="protein_similarity":
        settings = lock["protein_similarity"]
        result["alignment_audit"] = alignment_audit(training,test,Path(split)/"alignment_audit",settings["identity_threshold"],settings["coverage_threshold"])
        if result["alignment_audit"]["cross_threshold_matches"]:
            raise ValueError("Protein alignment leakage detected; split requires correction before interpreting results")
    write_json(out/"metrics.json",result)
    return result


def summarize(results,out):
    records = []
    for path in Path(results).rglob("metrics.json"):
        result = load_json(path)
        rank = result.get("within_target_ranking",{})
        records.append({"split":result["split_mode"],"model":result["model"],"seed":result["model_seed"],
                        **result["test_metrics"],"applicability_fraction":result["applicability_fraction"],
                        "rank_eligible_targets":rank.get("eligible_targets"),"within_target_macro_spearman":rank.get("macro_spearman_mean")})
    if not records:
        raise ValueError("No completed evaluation metrics")
    frame = pd.DataFrame(records).sort_values(["split","model","seed"])
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    frame.to_csv(out/"all_runs.csv",index=False)
    group = frame.groupby(["split","model"])[["rmse_pkd","mae_pkd","spearman","within_target_macro_spearman","applicability_fraction"]].agg(["mean","std"])
    group.to_csv(out/"summary.csv")
    write_json(out/"summary.json",{"runs":records,"interpretation":"Seed variability only; not a confidence interval. No external validation or food-specific efficacy claim.","completed_runs":len(records)})
    return records
