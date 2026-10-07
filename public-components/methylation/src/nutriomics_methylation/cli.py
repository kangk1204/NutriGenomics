from __future__ import annotations

import argparse
from pathlib import Path

from .io import fetch, prepare
from .model import train
from .evaluation import evaluate, export


def main(argv=None):
    parser = argparse.ArgumentParser(description="Research-only public hypertension DNAm pipeline")
    parser.add_argument("command", choices=["fetch", "prepare", "train", "evaluate", "export", "all", "transfer"])
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--cores", type=int, default=12)
    parser.add_argument("--spaces", nargs="+", choices=["common", "full"], default=["common", "full"])
    parser.add_argument("--methods", nargs="+", choices=["nested5", "loocv"], default=["nested5", "loocv"])
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    args = parser.parse_args(argv)
    root = args.root.resolve()
    commands = ["fetch", "prepare", "train", "evaluate", "export"] if args.command == "all" else [args.command]
    for command in commands:
        if command == "fetch":
            fetch(root)
        elif command == "prepare":
            prepare(root)
        elif command == "train":
            train(root, spaces=args.spaces, methods=args.methods, cores=args.cores, seed=args.seed)
        elif command == "evaluate":
            evaluate(root, args.bootstrap_replicates)
        elif command == "export":
            export(root)
        elif command == "transfer":
            from .transfer import reciprocal
            reciprocal(root, cores=args.cores, seed=args.seed, bootstraps=args.bootstrap_replicates)


if __name__ == "__main__":
    main()
