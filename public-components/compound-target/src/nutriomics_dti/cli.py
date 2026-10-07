import argparse
from pathlib import Path
import pandas as pd
from .io import write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description="Measured Kd regression with provenance-aware locked holdouts")
    commands = parser.add_subparsers(dest="command",required=True)
    fetch = commands.add_parser("fetch")
    fetch.add_argument("--out",required=True)
    fetch.add_argument("--subset",choices=["all","articles"],default="all")
    fetch.add_argument("--no-assays",action="store_true")
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--raw",required=True)
    prepare.add_argument("--assay-map")
    prepare.add_argument("--out",required=True)
    prepare.add_argument("--species",default=None)
    split = commands.add_parser("split")
    split.add_argument("--prepared",required=True)
    split.add_argument("--out",required=True)
    split.add_argument("--modes",nargs="+",default=["random","cold_compound","cold_target","scaffold","protein_similarity"])
    split.add_argument("--seeds",nargs="+",type=int,default=[17,29,43])
    split.add_argument("--protein-identity",type=float,default=0.4)
    split.add_argument("--protein-coverage",type=float,default=0.8)
    train = commands.add_parser("train")
    train.add_argument("--split",required=True)
    train.add_argument("--out",required=True)
    train.add_argument("--model",choices=["baseline","cnn"],default="baseline")
    train.add_argument("--seed",type=int,default=17)
    train.add_argument("--epochs",type=int,default=40)
    train.add_argument("--patience",type=int,default=6)
    train.add_argument("--trees",type=int,default=150)
    train.add_argument("--device",default="auto")
    evaluate = commands.add_parser("evaluate")
    evaluate.add_argument("--split",required=True)
    evaluate.add_argument("--model-dir",required=True)
    evaluate.add_argument("--out",required=True)
    evaluate.add_argument("--device",default="auto")
    predict = commands.add_parser("predict")
    predict.add_argument("--input",required=True,help="CSV with smiles and protein_sequence; no labels required")
    predict.add_argument("--model-dir",required=True)
    predict.add_argument("--out",required=True)
    predict.add_argument("--device",default="auto")
    summarize = commands.add_parser("summarize")
    summarize.add_argument("--results",required=True)
    summarize.add_argument("--out",required=True)
    arguments = parser.parse_args(argv)
    if arguments.command=="fetch":
        from .sources import fetch
        result = fetch(arguments.out,arguments.subset,not arguments.no_assays)
    elif arguments.command=="prepare":
        from .data import prepare
        result = prepare(arguments.raw,arguments.out,arguments.assay_map,arguments.species)
    elif arguments.command=="split":
        from .splits import lock_splits
        result = lock_splits(arguments.prepared,arguments.out,arguments.modes,arguments.seeds,arguments.protein_identity,arguments.protein_coverage)
    elif arguments.command=="train":
        from .models import train
        result = train(arguments.split,arguments.out,arguments.model,arguments.seed,arguments.epochs,arguments.patience,arguments.device,arguments.trees)
    elif arguments.command=="evaluate":
        from .models import evaluate
        result = evaluate(arguments.split,arguments.model_dir,arguments.out,arguments.device)
    elif arguments.command=="predict":
        from .models import predict
        from .data import structure,AA
        frame = pd.read_csv(arguments.input,keep_default_na=False)
        if not {"smiles","protein_sequence"}<=set(frame):
            parser.error("Input needs smiles and protein_sequence")
        frame["smiles"] = [structure(value)[0] for value in frame.smiles]
        frame["protein_sequence"] = frame.protein_sequence.str.upper()
        if any(len(seq)<20 or not set(seq)<=AA for seq in frame.protein_sequence):
            parser.error("Invalid protein sequence")
        output = predict(arguments.model_dir,frame,arguments.device)
        Path(arguments.out).parent.mkdir(parents=True,exist_ok=True)
        output.to_csv(arguments.out,index=False)
        result = {"predictions":len(output),"output":arguments.out,"status":"predictions, not measured evidence"}
    elif arguments.command=="summarize":
        from .models import summarize
        result = summarize(arguments.results,arguments.out)
    print(result,flush=True)


if __name__=="__main__":
    main()
