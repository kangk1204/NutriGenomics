"""Bounded baseline / DearDTI training using development-only file loaders."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from .features import features
from .models import metrics
from .io import sha256, utc_now, write_json, load_json
from .prospective_kd import development, consume_test, freeze, MODES


def code_hash():
    root = Path(__file__).parent
    inputs = [root/name for name in ("native_kd_runner.py","deardti_kd.py","prospective_kd.py","source_ledger_v2.py","vendor/deardti_encoders.py","vendor/deardti_fusion.py")]
    return {str(path.relative_to(root)):sha256(path) for path in inputs}


def train(split, out, model_name, seed=20261003, epochs=20, patience=4, device="cpu", threads=2, memory_fraction=.36):
    (training,validation),lock = development(split)
    if model_name not in lock["declared_models"]:
        raise ValueError("Model not declared before lock")
    out = Path(out)
    if (out/"training.json").exists():
        raise FileExistsError("Frozen training output already exists")
    out.mkdir(parents=True,exist_ok=True)
    started = time.monotonic()
    protocol = {"started_utc":utc_now(),"model":model_name,"split_lock_sha256":sha256(Path(split)/"lock.json"),"code_sha256":code_hash(),"fresh_weights":True,"uses_historical_weights":False,"training_loader":"train and validation only; test is never deserialized", "locked_test_sha256":lock["files"]["test"]["sha256"],"seed":seed,"device_requested":device,"threads":threads}
    # Record all bounded settings before the first fit, including interrupted fits.
    if model_name == "extratrees":
        config = {"trees":100,"max_depth":24,"min_samples_leaf_candidates":[2,4],"max_features":.7}
    elif model_name == "ridge":
        config = {"alpha_candidates":[1.,10.,100.],"scale":"training-only StandardScaler"}
    else:
        config = {"hidden":64,"gine_layers":3,"jk":True,"protein_encoder":"fresh AA-CNN 3layers", "protein_cap":1024,"protein_max_tokens":128,"max_atoms":128,"batch_size":32,"epochs_max":epochs,"patience":patience,"lr":.001,"weight_decay":.0001,"gradient_clip":5.,"cuda_memory_fraction":memory_fraction}
    protocol["config"] = config
    if (out/"training_started.json").exists():
        raise FileExistsError("Interrupted run retained; use a new directory")
    write_json(out/"training_started.json",protocol)
    history = []
    runtime = {"python":platform.python_version(),"platform":platform.platform(),"numpy":np.__version__}
    if model_name in {"extratrees","ridge"}:
        import sklearn
        runtime["sklearn"] = sklearn.__version__
        xtrain,xval = features(training),features(validation)
        candidates = config.get("min_samples_leaf_candidates",config.get("alpha_candidates"))
        best_score = float("inf")
        best_model = None
        for value in candidates:
            if model_name == "extratrees":
                model = ExtraTreesRegressor(n_estimators=100,max_depth=24,min_samples_leaf=value,max_features=.7,random_state=seed,n_jobs=threads)
            else:
                model = make_pipeline(StandardScaler(),Ridge(alpha=value))
            model.fit(xtrain,training.pkd)
            predicted = model.predict(xval)
            score = metrics(validation.pkd,predicted)
            history.append({"parameter":value,"validation":score})
            if score["rmse_pkd"] < best_score:
                best_score,best_model,selected_parameter = score["rmse_pkd"],model,value
        joblib.dump(best_model,out/"model.joblib")
        prediction = best_model.predict(xval)
        config["selected_parameter"] = selected_parameter
    elif model_name == "deardti_gine_aacnn":
        import random
        import torch
        from .deardti_kd import NativeKdGINE, PairDataset, collate, predict
        random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
        torch.set_num_threads(threads)
        torch.use_deterministic_algorithms(True,warn_only=True)
        torch.backends.cudnn.benchmark = False
        if device == "cuda":
            if not torch.cuda.is_available():
                raise RuntimeError("CUDA requested but unavailable")
            torch.cuda.manual_seed_all(seed)
            torch.cuda.set_per_process_memory_fraction(memory_fraction)
        model = NativeKdGINE().to(device)
        mean,scale = float(training.pkd.mean()),max(float(training.pkd.std()),1e-6)
        loader = torch.utils.data.DataLoader(PairDataset(training,mean,scale),batch_size=32,shuffle=True,collate_fn=collate,generator=torch.Generator().manual_seed(seed),num_workers=0)
        optimizer = torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.0001)
        best_score,best_state,stale = float("inf"),None,0
        for epoch in range(1,epochs+1):
            model.train(); losses=[]
            for batch in loader:
                batch = {key:value.to(device) if isinstance(value,torch.Tensor) else value for key,value in batch.items()}
                optimizer.zero_grad(set_to_none=True)
                prediction = model(batch)
                loss = torch.nn.functional.mse_loss(prediction,batch["labels"])
                loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),5.)
                optimizer.step(); losses.append(float(loss.detach().cpu()))
            prediction = predict(model,validation,device,mean,scale)
            score = metrics(validation.pkd,prediction)
            history.append({"epoch":epoch,"train_normalized_mse":float(np.mean(losses)),"validation":score})
            write_json(out/"progress.json",{"epoch":epoch,"history":history,"elapsed_seconds":time.monotonic()-started})
            print(f"{lock['mode']} epoch={epoch} valRMSE={score['rmse_pkd']:.6f}",flush=True)
            if score["rmse_pkd"] < best_score:
                best_score,best_state,stale = score["rmse_pkd"],{key:value.detach().cpu().clone() for key,value in model.state_dict().items()},0
                config["selected_epoch"] = epoch
            else:
                stale += 1
                if stale >= patience:
                    break
        if best_state is None:
            raise ValueError("No valid validation checkpoint")
        model.load_state_dict(best_state)
        torch.save({"state_dict":best_state,"mean":mean,"scale":scale},out/"deardti_kd.pt")
        prediction = predict(model,validation,device,mean,scale)
        runtime.update({"torch":torch.__version__,"cuda_build":torch.version.cuda,"gpu":torch.cuda.get_device_name(0) if device=="cuda" else None,"peak_cuda_allocated_bytes":torch.cuda.max_memory_allocated() if device=="cuda" else None})
        config["normalization"] = {"training_mean":mean,"training_std":scale}
    else:
        raise ValueError(model_name)
    training[["compound_id","target_id","smiles","protein_sequence","pair_id"]].drop_duplicates().to_csv(out/"training_entities.csv.gz",index=False)
    record = {**protocol,"completed_utc":utc_now(),"elapsed_seconds":time.monotonic()-started,"config":config,"runtime":runtime,"history":history,"validation_metrics":metrics(validation.pkd,prediction),"model_sha256":sha256(out/("deardti_kd.pt" if model_name=="deardti_gine_aacnn" else "model.joblib")),"selection":"development validation only; test attempt remains unconsumed"}
    write_json(out/"training.json",record)
    return record


def predict_model(model_dir,frame,device="cpu"):
    model_dir = Path(model_dir)
    record = load_json(model_dir/"training.json")
    model = record["model"]
    path = model_dir/("deardti_kd.pt" if model=="deardti_gine_aacnn" else "model.joblib")
    if sha256(path)!=record["model_sha256"]:
        raise ValueError("Frozen model changed")
    if model != "deardti_gine_aacnn":
        return joblib.load(path).predict(features(frame))
    import torch
    from .deardti_kd import NativeKdGINE,predict
    torch.set_num_threads(record["threads"])
    if device=="cuda":
        torch.cuda.set_per_process_memory_fraction(record["config"]["cuda_memory_fraction"])
    checkpoint = torch.load(path,map_location="cpu",weights_only=True)
    model = NativeKdGINE().to(device)
    model.load_state_dict(checkpoint["state_dict"])
    return predict(model,frame,device,checkpoint["mean"],checkpoint["scale"])


def evaluate(split,model_dir,out,device="cpu"):
    model_dir,out = Path(model_dir),Path(out)
    record = load_json(model_dir/"training.json")
    if record["split_lock_sha256"] != sha256(Path(split)/"lock.json"):
        raise ValueError("Different locked split")
    # Consume the attempt before deserializing any test labels.
    test,lock = consume_test(split,record["model"],out)
    predictions = predict_model(model_dir,test,device)
    result = test[["measurement_id","pair_id","compound_id","target_id","source_group","pkd"]].copy()
    result["predicted_pkd"] = predictions
    result["evidence_status"] = "prediction; not measured interaction or food efficacy"
    result.to_csv(out/"test_predictions.csv.gz",index=False)
    per_source = []
    for source,subset in result.groupby("source_group"):
        per_source.append({"source_group":source,**metrics(subset.pkd,subset.predicted_pkd)})
    pd.DataFrame(per_source).to_csv(out/"source_metrics.csv",index=False)
    foodmin = {"n":30,"compounds":5,"targets":3,"source_groups":3}
    report = {"evaluated_utc":utc_now(),"model":record["model"],"mode":lock["mode"],"seed":lock["seed"],"protocol_id":lock["protocol_id"],"test_sha256":lock["files"]["test"]["sha256"],"model_sha256":record["model_sha256"],"metrics":metrics(test.pkd,predictions),"source_groups":len(per_source),"overlap_audit":lock["overlap_counts"],"evaluation_type":lock["evaluation_type"],"historical_external_independence":False,"food_independent_evaluation":"unevaluated until complete historical ledger and minimum food cohort guards pass","food_minimum":foodmin,"test_used_for_model_selection":False}
    write_json(out/"metrics.json",report)
    write_json(Path(split)/f"test_result_{record['model']}.json",{"completed_utc":utc_now(),"metrics_sha256":sha256(out/"metrics.json"),"predictions_sha256":sha256(out/"test_predictions.csv.gz")})
    return report


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command",required=True)
    ledger = sub.add_parser("ledger")
    ledger.add_argument("--prepared",required=True); ledger.add_argument("--prepare-audit"); ledger.add_argument("--out",required=True)
    lock = sub.add_parser("freeze")
    lock.add_argument("--prepared",required=True); lock.add_argument("--out",required=True); lock.add_argument("--modes",nargs="+",default=list(MODES)); lock.add_argument("--threads",type=int,default=2); lock.add_argument("--components")
    fit = sub.add_parser("train")
    fit.add_argument("--split",required=True); fit.add_argument("--out",required=True); fit.add_argument("--model",dest="model_name",choices=("extratrees","ridge","deardti_gine_aacnn"),required=True); fit.add_argument("--epochs",type=int,default=20); fit.add_argument("--patience",type=int,default=4); fit.add_argument("--device",default="cpu"); fit.add_argument("--threads",type=int,default=2)
    assess = sub.add_parser("evaluate")
    assess.add_argument("--split",required=True); assess.add_argument("--model-dir",required=True); assess.add_argument("--out",required=True); assess.add_argument("--device",default="cpu")
    args = vars(parser.parse_args()); command = args.pop("command")
    if command=="ledger":
        from .source_ledger_v2 import build
        print(build(**args))
    else:
        print({"freeze":freeze,"train":train,"evaluate":evaluate}[command](**args))


if __name__=="__main__":
    main()
