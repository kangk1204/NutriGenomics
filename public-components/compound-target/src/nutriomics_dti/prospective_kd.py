"""New, source-closed native Kd benchmark, distinct from all old test recipes.

Training loaders never deserialize the test file. Each declared model gets one
final-test attempt (including failures); test scores cannot choose a model.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from .io import sha256, utc_now, write_json, load_json
from .source_ledger_v2 import native_guards, close_aliases
from .splits import assert_no_overlap

MODES = ("random", "cold_compound", "cold_target", "scaffold", "protein_similarity", "source")
MODELS = ("extratrees", "ridge", "deardti_gine_aacnn")
SEED = 20261003


def source_overlap(frames):
    sets = [{alias for text in frame.source_aliases_json for alias in json.loads(text)} for frame in frames]
    overlaps = [len(sets[a]&sets[b]) for a, b in ((0,1), (0,2), (1,2))]
    if any(overlaps):
        raise ValueError("Canonical source alias leaks across partitions")
    return overlaps


def partition(frame, mode, seed, min_rows=20):
    keys = {"random": "pair_id", "cold_compound": "compound_id", "cold_target": "target_id", "scaffold": "scaffold", "protein_similarity": "protein_cluster", "source": "source_group"}
    key = keys[mode]
    groups = np.asarray(sorted(frame[key].unique()))
    if len(groups) < 5:
        raise ValueError("Too few groups for sealed partitions")
    shuffled = np.random.default_rng(seed).permutation(groups)
    ntest = max(1, int(len(groups)*.15))
    nval = max(1, int(len(groups)*.15))
    test = frame[frame[key].isin(shuffled[:ntest])].copy()
    validation = frame[frame[key].isin(shuffled[ntest:ntest+nval])].copy()
    train = frame[~frame[key].isin(shuffled[:ntest+nval])].copy()
    badval = validation.source_group.isin(test.source_group) | validation.pair_id.isin(test.pair_id)
    qval = validation[badval].assign(reason="test_source_or_pair")
    validation = validation[~badval]
    heldsource = set(test.source_group)|set(validation.source_group)
    heldpair = set(test.pair_id)|set(validation.pair_id)
    badtrain = train.source_group.isin(heldsource)|train.pair_id.isin(heldpair)
    qtrain = train[badtrain].assign(reason="heldout_source_or_pair")
    train = train[~badtrain]
    frames = (train, validation, test)
    overlaps = assert_no_overlap(*frames, mode)
    overlaps["canonical_source_alias"] = source_overlap(frames)
    if min(map(len, frames)) < min_rows:
        raise ValueError(f"Source-purged partitions too small: {[len(f) for f in frames]}")
    return dict(zip(("train", "validation", "test"), frames)), pd.concat((qval, qtrain)), overlaps


def freeze(prepared, out, modes=MODES, seeds=(SEED,), threads=4, components=None):
    from .protein_split_v2 import search, edge_components
    out = Path(out)
    if (out/"protocol.json").exists():
        raise FileExistsError("Protocol already frozen; do not regenerate test membership")
    frame = pd.read_csv(prepared, keep_default_na=False)
    native_guards(frame)
    frame = close_aliases(frame)
    # Dense graph encoder never silently drops atoms. Apply same eligibility to controls.
    atoms = frame.smiles.map(lambda value: Chem.MolFromSmiles(value).GetNumAtoms())
    excluded = frame[atoms>128].copy()
    frame = frame[atoms<=128].copy()
    out.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(out/"encoder_excluded.csv.gz", index=False)
    sequence_settings = {}
    if "protein_similarity" in modes:
        proteins = frame[["target_id", "protein_sequence"]].drop_duplicates()
        if components:
            comp = pd.read_csv(components, sep="\t", keep_default_na=False)
            mapping = dict(zip(comp.target_id, comp.protein_cluster))
            if not set(proteins.target_id)<=set(mapping):
                raise ValueError("Existing components do not cover every current sequence")
            sequence_settings = {"components_input_sha256": sha256(components), "identity": .4, "coverage": .8}
        else:
            edges, provenance = search(proteins, proteins, out/"sequence_all_vs_all", .4, .8, threads)
            mapping = edge_components(proteins.target_id.tolist(), edges)
            sequence_settings = {"all_vs_all": provenance, "identity": .4, "coverage": .8}
        frame["protein_cluster"] = frame.target_id.map(mapping)
        pd.DataFrame(sorted(mapping.items()), columns=["target_id", "protein_cluster"]).to_csv(out/"protein_components.tsv", sep="\t", index=False)
    else:
        frame["protein_cluster"] = frame.target_id
    summary = []
    for mode in modes:
        if mode not in MODES:
            raise ValueError(mode)
        for seed in seeds:
            directory = out/mode/str(seed)
            if (directory/"lock.json").exists():
                raise FileExistsError("Split lock already exists")
            partitions, quarantine, overlaps = partition(frame, mode, seed)
            gates = {}
            if mode == "protein_similarity":
                for a,b in (("train","validation"),("validation","train"),("train","test"),("test","train"),("validation","test"),("test","validation")):
                    found, record = search(partitions[a], partitions[b], out/"sequence_gates"/str(seed)/f"{a}_to_{b}", .4, .8, threads)
                    if found:
                        raise ValueError("Cross-threshold sequence alignment: do not release this split")
                    gates[f"{a}_to_{b}"] = record
            directory.mkdir(parents=True, exist_ok=True)
            files = {}
            for name, subset in partitions.items():
                path = directory/f"{name}.csv.gz"
                subset.sort_values("measurement_id").to_csv(path, index=False)
                files[name] = {"rows": len(subset), "sha256": sha256(path), "pairs": int(subset.pair_id.nunique()), "source_groups": int(subset.source_group.nunique())}
            quarantine.to_csv(directory/"quarantine.csv.gz", index=False)
            lock = {"protocol_id": "native_kd_gine_20261003_v1", "frozen_utc": utc_now(), "mode": mode, "seed": seed, "prepared_sha256": sha256(prepared), "files": files, "overlap_counts": overlaps, "sequence_alignment_gates": gates, "declared_models": list(MODELS), "historical_external_independence": False, "evaluation_type": "new internal source-purged split; earlier corpus exposure exists; not independent validation of historical models", "partition": "70/15/15 groups before source/pair quarantine", "quarantined": len(quarantine), "molecule_max_atoms": 128}
            write_json(directory/"lock.json", lock)
            summary.append({"mode": mode, "seed": seed, **{name: val["rows"] for name,val in files.items()}, "quarantined":len(quarantine)})
            print(summary[-1], flush=True)
    protocol = {"protocol_id": "native_kd_gine_20261003_v1", "frozen_utc": utc_now(), "input_sha256": sha256(prepared), "eligible_records": len(frame), "excluded_gt128_atoms": len(excluded), "endpoint": "exact positive native Kd nM; pKd=-log10(Kd[M])", "models": list(MODELS), "seeds": list(seeds), "modes": list(modes), "summary": summary, "protein_similarity": sequence_settings, "selection": "Only validation selects baseline parameters or neural checkpoint; report all declared test models", "fresh_supervised_weights": True, "sequence_encoder": "fresh amino acid embedding + 3 layer CNN; no pretrained embeddings", "historical_independence": "Not established without immutable complete historical ledger", "test_attempt_rule": "One test attempt per declared model/seed/split; any failed attempt is retained and consumes the attempt"}
    write_json(out/"protocol.json", protocol)
    return protocol


def development(directory):
    directory = Path(directory)
    lock = load_json(directory/"lock.json")
    frames = []
    for name in ("train", "validation"):
        path = directory/f"{name}.csv.gz"
        if sha256(path) != lock["files"][name]["sha256"]:
            raise ValueError("Development input differs from frozen lock")
        frames.append(pd.read_csv(path, keep_default_na=False))
    return frames, lock


def consume_test(directory, model, out):
    directory, out = Path(directory), Path(out)
    lock = load_json(directory/"lock.json")
    if model not in lock["declared_models"]:
        raise ValueError("Model was not declared before test freezing")
    out.mkdir(parents=True, exist_ok=True)
    attempt = directory/f"test_attempt_{model}.json"
    record = {"started_utc": utc_now(), "model": model, "lock_sha256": sha256(directory/"lock.json"), "test_sha256": lock["files"]["test"]["sha256"], "out": str(out), "status": "attempt_consumed_before_label_access"}
    # Exclusive creation protects concurrent readers and retries after failure.
    with attempt.open("x", encoding="utf-8") as stream:
        json.dump(record, stream, indent=2)
    test = directory/"test.csv.gz"
    if sha256(test) != record["test_sha256"]:
        raise ValueError("Sealed test input changed; attempt consumed")
    return pd.read_csv(test, keep_default_na=False), lock
