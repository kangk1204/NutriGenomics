"""Reported-edge connected components and pretraining cross-fold alignment gates.

Version 1 used MMseqs easy-cluster, which is a cascaded clustering workflow and
does not establish that every qualifying search edge stays in one component.
Version 2 unions all threshold-qualified all-vs-all search edges, then checks
each direction of every partition pair before releasing any new training lock.
Old files/results are preserved. Partition corrections never consult Kd labels.
"""
import argparse
import shutil
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from .io import load_json,sha256,utc_now,write_json
from .splits import UnionFind,write_fasta,assert_no_overlap


def qualifying_edges(path,identity=.4,coverage=.8):
    path = Path(path)
    if not path.stat().st_size:
        return []
    frame = pd.read_csv(path,sep="\t",names=["query","target","identity","query_coverage","target_coverage"])
    numeric = frame[["identity","query_coverage","target_coverage"]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all() or (numeric<0).any() or (numeric>1+1e-6).any():
        raise ValueError("MMseqs fident/qcov/tcov must be finite fractions in [0,1]")
    keep = frame[(frame.identity>=identity-1e-6)&(frame.query_coverage>=coverage-1e-6)&(frame.target_coverage>=coverage-1e-6)]
    return sorted(set(zip(keep["query"],keep["target"])))


def edge_components(identifiers,edges):
    groups = UnionFind(identifiers)
    known = set(identifiers)
    for query,target in edges:
        if query not in known or target not in known:
            raise ValueError("Search edge references a target outside the frozen cohort")
        groups.union(query,target)
    return {identifier:groups.find(identifier) for identifier in identifiers}


def search(query,reference,out,identity=.4,coverage=.8,threads=6):
    out = Path(out)
    out.mkdir(parents=True,exist_ok=True)
    result = out/"alignments.tsv"
    query_fasta,reference_fasta = out/"query.fasta",out/"reference.fasta"
    write_fasta(query,query_fasta)
    write_fasta(reference,reference_fasta)
    max_sequences = int(reference.target_id.nunique())
    if max_sequences<1 or query.target_id.nunique()<1:
        raise ValueError("Cannot audit an empty target partition")
    # Exhaustive-search avoids representative-only and top-300 candidate caps.
    # Explicit post-filtering uses the same fident and bidirectional coverage.
    command = ["mmseqs","easy-search",str(query_fasta),str(reference_fasta),str(result),str(out/"tmp"),
               "--min-seq-id",str(identity),"-c",str(coverage),"--cov-mode","0","--alignment-mode","3",
               "--exhaustive-search","1","--max-seqs",str(max_sequences),"-e","1000000",
               "-s","7.5","--threads",str(threads),"--format-output","query,target,fident,qcov,tcov"]
    with open(out/"search.log","w") as stream:
        subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True)
    edges = qualifying_edges(result,identity,coverage)
    record = {"command":command,"query_fasta_sha256":sha256(query_fasta),"reference_fasta_sha256":sha256(reference_fasta),
              "alignments_sha256":sha256(result),"qualified_edges":len(edges),"max_sequences":max_sequences,
              "identity_threshold":identity,"bidirectional_coverage_threshold":coverage,
              "interpretation":"MMseqs exhaustive candidate search and its reported local alignments; not a proof over every possible global alignment"}
    write_json(out/"provenance.json",record)
    return edges,record


def partition(frame,mapping,seed):
    frame = frame.copy()
    frame["protein_cluster"] = frame.target_id.map(mapping)
    if frame.protein_cluster.isna().any():
        raise ValueError("Missing target in v2 components")
    groups = sorted(frame.protein_cluster.unique())
    if len(groups)<5:
        raise ValueError("Too few connected components for three valid partitions")
    shuffled = np.random.default_rng(seed).permutation(groups)
    ntest,nval = max(1,int(len(groups)*.15)),max(1,int(len(groups)*.15))
    test_groups,val_groups = set(shuffled[:ntest]),set(shuffled[ntest:ntest+nval])
    test = frame[frame.protein_cluster.isin(test_groups)].copy()
    validation = frame[frame.protein_cluster.isin(val_groups)].copy()
    train = frame[~frame.protein_cluster.isin(test_groups|val_groups)].copy()
    excluded_validation = validation[validation.source_group.isin(test.source_group)|validation.pair_id.isin(test.pair_id)]
    validation = validation.drop(excluded_validation.index)
    held_sources,held_pairs = set(test.source_group)|set(validation.source_group),set(test.pair_id)|set(validation.pair_id)
    excluded_train = train[train.source_group.isin(held_sources)|train.pair_id.isin(held_pairs)]
    train = train.drop(excluded_train.index)
    audit = assert_no_overlap(train,validation,test,"protein_similarity")
    if min(len(train),len(validation),len(test))<20:
        raise ValueError("V2 provenance-purged split has fewer than 20 measurements in a partition")
    quarantine = pd.concat([excluded_validation.assign(reason="test_source_or_pair"),excluded_train.assign(reason="heldout_source_or_pair")])
    return {"train":train,"validation":validation,"test":test},quarantine,audit


def prepare(prepared,out,previous=None,seeds=(17,29,43),identity=.4,coverage=.8,threads=6,max_rounds=5):
    out = Path(out)
    if (out/"split_summary.json").exists():
        raise FileExistsError("V2 splits already frozen; retain them and choose a new output directory")
    if not shutil.which("mmseqs"):
        raise RuntimeError("MMseqs2 is required")
    out.mkdir(parents=True,exist_ok=True)
    frame = pd.read_csv(prepared,keep_default_na=False)
    proteins = frame[["target_id","protein_sequence"]].drop_duplicates().sort_values("target_id")
    if proteins.target_id.duplicated().any():
        raise ValueError("Target ID refers to conflicting protein sequences")
    identifiers = proteins.target_id.tolist()
    edges,initial_search = search(proteins,proteins,out/"all_vs_all",identity,coverage,threads)
    all_edges = set(edges)
    diagnostics = []
    final = None
    for round_number in range(1,max_rounds+1):
        mapping = edge_components(identifiers,sorted(all_edges))
        candidates,new_edges = {},set()
        for seed in seeds:
            partitions,quarantine,overlaps = partition(frame,mapping,seed)
            gates = {}
            for a,b in (("train","validation"),("validation","train"),("train","test"),("test","train"),("validation","test"),("test","validation")):
                found,record = search(partitions[a],partitions[b],out/"gates"/f"round_{round_number}"/str(seed)/f"{a}_to_{b}",identity,coverage,threads)
                gates[f"{a}_to_{b}"] = record
                new_edges.update(found)
            candidates[seed] = (partitions,quarantine,overlaps,gates)
        new_edges -= all_edges
        diagnostics.append({"round":round_number,"components":len(set(mapping.values())),"additional_cross_edges":len(new_edges)})
        write_json(out/"refinement_diagnostics.json",diagnostics)
        if not new_edges:
            if any(gate["qualified_edges"] for _,_,_,gates in candidates.values() for gate in gates.values()):
                raise ValueError("Cross-fold alignment remains despite component assignment; no locks released")
            final = (mapping,candidates)
            break
        all_edges.update(new_edges)
    if final is None:
        raise ValueError("Cross-fold refinement did not converge; no v2 locks released")
    mapping,candidates = final
    pd.DataFrame(sorted(mapping.items()),columns=["target_id","protein_cluster"]).to_csv(out/"protein_components.tsv",sep="\t",index=False)
    version = subprocess.check_output(["mmseqs","version"],text=True).strip()
    settings = {"protocol_version":"v2_reported_edge_connected_components","metric":"MMseqs fident plus bidirectional coverage",
                "identity_threshold":identity,"coverage_threshold":coverage,"mmseqs_version":version,
                "all_vs_all":initial_search,"refinement":diagnostics,"components":len(set(mapping.values())),
                "components_sha256":sha256(out/"protein_components.tsv"),"qualified_all_edges":len(all_edges),
                "selection":"Protein sequences/topology only; no affinity-label selection"}
    summaries = []
    invalid = []
    for seed,(partitions,quarantine,overlaps,gates) in candidates.items():
        directory = out/"protein_similarity"/str(seed)
        directory.mkdir(parents=True,exist_ok=True)
        if (directory/"lock.json").exists():
            raise FileExistsError("Cannot overwrite an existing v2 lock")
        files = {}
        for name,subset in partitions.items():
            path = directory/f"{name}.csv.gz"
            subset.sort_values("measurement_id").to_csv(path,index=False)
            files[name] = {"rows":len(subset),"sha256":sha256(path),"pairs":subset.pair_id.nunique(),"source_groups":subset.source_group.nunique()}
        quarantine.to_csv(directory/"quarantine.csv.gz",index=False)
        lock = {"locked_utc":utc_now(),"mode":"protein_similarity","seed":seed,"protocol_version":settings["protocol_version"],
                "prepared_sha256":sha256(prepared),"files":files,"overlap_counts":overlaps,"protein_similarity":settings,
                "pretraining_alignment_gates":gates,"quarantined_train":int((quarantine.reason=="heldout_source_or_pair").sum()),
                "quarantined_validation":int((quarantine.reason=="test_source_or_pair").sum()),
                "raw_partition_ratio":"70/15/15 connected components before provenance quarantine",
                "evaluation_type":"corrected internal holdout; v1 topology correction after a failed audit, not independent external validation"}
        write_json(directory/"lock.json",lock)
        summaries.append({"mode":"protein_similarity","seed":seed,"protocol_version":settings["protocol_version"],**{name:row["rows"] for name,row in files.items()},
                          "all_six_direction_cross_threshold_matches":0})
        if previous:
            old = Path(previous)/"protein_similarity"/str(seed)/"lock.json"
            if old.exists():
                invalid.append({"path":str(old),"sha256":sha256(old),"status":"invalid_v1_protein_similarity_protocol",
                                "reason":"At least one qualifying cross-fold alignment was observed at v1 seed17; all v1 protein-similarity locks superseded, files preserved"})
    write_json(out/"invalid_v1_protocol.json",{"recorded_utc":utc_now(),"old_locks":invalid,"raw_and_other_four_split_modes_unchanged":True})
    write_json(out/"protocol.json",settings)
    write_json(out/"split_summary.json",summaries)
    return summaries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepared",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--previous",type=Path)
    parser.add_argument("--seeds",type=int,nargs="+",default=[17,29,43])
    parser.add_argument("--threads",type=int,default=6)
    args = parser.parse_args()
    for row in prepare(args.prepared,args.out,args.previous,tuple(args.seeds),threads=args.threads):
        print(row,flush=True)


if __name__=="__main__":
    main()
