"""Locked, provenance-purged splits. Structural similarity is reported, not hidden."""
import hashlib
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from .io import sha256, utc_now, write_json

MODES = ("random", "cold_compound", "cold_target", "scaffold", "protein_similarity")


class UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}
    def find(self, value):
        while value != self.parent[value]:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value
    def union(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def kmers(sequence, k=3):
    return {sequence[index:index+k] for index in range(len(sequence)-k+1)}


def kmer_jaccard(a, b):
    return len(a & b) / len(a | b) if a or b else 1.0


def protein_clusters(frame, threshold=0.5):
    proteins = frame.drop_duplicates("target_id").sort_values("target_id")
    identifiers = proteins.target_id.tolist()
    sets = [kmers(sequence) for sequence in proteins.protein_sequence]
    groups = UnionFind(identifiers)
    edges = 0
    for i, current in enumerate(sets):
        for j in range(i):
            # This upper bound avoids expensive impossible comparisons.
            other = sets[j]
            if min(len(current), len(other)) / max(len(current), len(other)) < threshold:
                continue
            if kmer_jaccard(current, other) >= threshold:
                groups.union(identifiers[i], identifiers[j])
                edges += 1
    return {identifier: groups.find(identifier) for identifier in identifiers}, edges


def write_fasta(frame,path):
    with open(path,"w",encoding="ascii") as stream:
        for row in frame.drop_duplicates("target_id").sort_values("target_id").itertuples():
            stream.write(f">{row.target_id}\n{row.protein_sequence}\n")


def mmseqs_clusters(frame,out,identity=0.4,coverage=0.8):
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    fasta = out/"proteins.fasta"
    write_fasta(frame,fasta)
    command = ["mmseqs","easy-cluster",str(fasta),str(out/"clusters"),str(out/"tmp"),
               "--min-seq-id",str(identity),"-c",str(coverage),"--cov-mode","0","--cluster-mode","1",
               "--alignment-mode","3","-s","7.5","--threads","6"]
    with open(out/"clustering.log","w") as stream:
        subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
    clusters = pd.read_csv(out/"clusters_cluster.tsv",sep="\t",names=["cluster","target"])
    mapping = dict(zip(clusters.target,clusters.cluster))
    if set(mapping)!=set(frame.target_id):
        raise ValueError("MMseqs did not cluster every target")
    version = subprocess.check_output(["mmseqs","version"],text=True).strip()
    write_json(out/"provenance.json",{"command":command,"version":version,"fasta_sha256":sha256(fasta)})
    return mapping,{"metric":"MMseqs2 aligned sequence identity; bidirectional coverage; connected-component clustering",
                    "identity_threshold":identity,"coverage_threshold":coverage,"mmseqs_version":version,"command":command}


def alignment_audit(train,test,out,identity=0.4,coverage=0.8):
    """Search heldout proteins against train at the prespecified alignment threshold."""
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    result_path = out/"cross_alignments.tsv"
    if not result_path.exists():
        write_fasta(train,out/"train.fasta")
        write_fasta(test,out/"test.fasta")
        command = ["mmseqs","easy-search",str(out/"test.fasta"),str(out/"train.fasta"),str(result_path),str(out/"tmp"),
                   "--min-seq-id",str(identity),"-c",str(coverage),"--cov-mode","0","--alignment-mode","3",
                   "-s","7.5","--threads","6","--format-output","query,target,fident,qcov,tcov"]
        with open(out/"search.log","w") as stream:
            subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
        write_json(out/"search_provenance.json",{"command":command,"identity_threshold":identity,"coverage_threshold":coverage})
    if result_path.stat().st_size:
        matches = pd.read_csv(result_path,sep="\t",names=["query","target","identity","query_coverage","target_coverage"])
        violating = matches[(matches.identity>=identity-1e-6)&(matches.query_coverage>=coverage-1e-6)&(matches.target_coverage>=coverage-1e-6)]
    else:
        violating = []
    return {"cross_threshold_matches":len(violating),"identity_threshold":identity,"bidirectional_coverage_threshold":coverage,
            "interpretation":"MMseqs2 search audit at sensitivity 7.5; not an exhaustive global-alignment proof"}


def assert_no_overlap(train, validation, test, mode):
    fields = ["pair_id", "source_group"]
    if mode == "cold_compound":
        fields.append("compound_id")
    if mode in ("cold_target", "protein_similarity"):
        fields.append("target_id")
    if mode == "scaffold":
        fields.append("scaffold")
    if mode == "protein_similarity":
        fields.append("protein_cluster")
    frames = (train, validation, test)
    report = {}
    for field in fields:
        sets = [set(frame[field]) for frame in frames]
        overlaps = [len(sets[i] & sets[j]) for i, j in ((0,1),(0,2),(1,2))]
        report[field] = overlaps
        if any(overlaps):
            raise ValueError(f"Leakage in {field}: train/val, train/test, val/test={overlaps}")
    return report


def lock_splits(prepared, out, modes=MODES, seeds=(17,29,43), protein_threshold=0.4, protein_coverage=0.8):
    out = Path(out)
    frame = pd.read_csv(prepared, keep_default_na=False)
    similarity = {"metric":"not requested"}
    cluster_map = {}
    if "protein_similarity" in modes:
        if not shutil.which("mmseqs"):
            raise RuntimeError("MMseqs2 is required for the primary protein similarity split; do not substitute k-mer similarity silently")
        cluster_map,similarity = mmseqs_clusters(frame,out/"protein_clustering",protein_threshold,protein_coverage)
    frame["protein_cluster"] = frame.target_id.map(cluster_map).fillna(frame.target_id)
    key_fields = {"random":"pair_id", "cold_compound":"compound_id", "cold_target":"target_id",
                  "scaffold":"scaffold", "protein_similarity":"protein_cluster"}
    summaries = []
    for mode in modes:
        key = key_fields[mode]
        groups = sorted(frame[key].unique())
        if len(groups) < 5:
            raise ValueError(f"Too few {mode} groups for train/validation/test: {len(groups)}")
        for seed in seeds:
            directory = out / mode / str(seed)
            if (directory / "lock.json").exists():
                raise FileExistsError(f"Split already locked: {directory}; use a new output directory")
            directory.mkdir(parents=True, exist_ok=True)
            permutation = np.random.default_rng(seed).permutation(groups)
            ntest = max(1, int(len(groups)*0.15))
            nval = max(1, int(len(groups)*0.15))
            test_groups = set(permutation[:ntest])
            val_groups = set(permutation[ntest:ntest+nval])
            test = frame[frame[key].isin(test_groups)].copy()
            validation = frame[frame[key].isin(val_groups)].copy()
            train = frame[~frame[key].isin(test_groups | val_groups)].copy()
            # Keep test fixed; purge source and pair relatives from development.
            excluded_val = validation[validation.source_group.isin(test.source_group) | validation.pair_id.isin(test.pair_id)]
            validation = validation.drop(excluded_val.index)
            held_sources = set(test.source_group) | set(validation.source_group)
            held_pairs = set(test.pair_id) | set(validation.pair_id)
            excluded_train = train[train.source_group.isin(held_sources) | train.pair_id.isin(held_pairs)]
            train = train.drop(excluded_train.index)
            overlaps = assert_no_overlap(train, validation, test, mode)
            if min(len(train), len(validation), len(test)) < 20:
                raise ValueError(f"Source-purged {mode}/{seed} has too few rows: {len(train)}, {len(validation)}, {len(test)}")
            files = {}
            for name, subset in (("train",train),("validation",validation),("test",test)):
                subset.sort_values("measurement_id").to_csv(directory / f"{name}.csv.gz", index=False)
                files[name] = {"rows":len(subset), "sha256":sha256(directory/f"{name}.csv.gz"),
                               "pairs":subset.pair_id.nunique(), "source_groups":subset.source_group.nunique()}
            pd.concat([excluded_val.assign(reason="test_source_or_pair"), excluded_train.assign(reason="heldout_source_or_pair")]).to_csv(directory/"quarantine.csv.gz", index=False)
            lock = {"locked_utc":utc_now(), "mode":mode, "seed":seed, "prepared_sha256":sha256(prepared),
                    "files":files, "overlap_counts":overlaps, "raw_partition_ratio":"70/15/15 groups before provenance quarantine",
                    "quarantined_train":len(excluded_train), "quarantined_validation":len(excluded_val),
                    "protein_similarity":similarity,
                    "evaluation_type":"internal locked holdout, not independent external validation"}
            write_json(directory/"lock.json",lock)
            summaries.append({"mode":mode,"seed":seed,**{name: entry["rows"] for name,entry in files.items()},
                              "quarantined":len(excluded_val)+len(excluded_train)})
            print(summaries[-1],flush=True)
    write_json(out/"split_summary.json",summaries)
    return summaries


def load_locked(directory):
    from .io import load_json
    directory = Path(directory)
    lock = load_json(directory/"lock.json")
    for name, entry in lock["files"].items():
        if sha256(directory/f"{name}.csv.gz") != entry["sha256"]:
            raise ValueError(f"Locked {name} file changed")
    frames = [pd.read_csv(directory/f"{name}.csv.gz",keep_default_na=False) for name in ("train","validation","test")]
    assert_no_overlap(*frames,lock["mode"])
    return frames, lock
