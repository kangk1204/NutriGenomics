from functools import lru_cache
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from .splits import kmers, kmer_jaccard

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024,includeChirality=True)


@lru_cache(maxsize=100000)
def fingerprint(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES")
    return MORGAN.GetFingerprint(molecule)


@lru_cache(maxsize=100000)
def drug_features(smiles):
    array = np.zeros(1024,dtype=np.float32)
    DataStructs.ConvertToNumpyArray(fingerprint(smiles),array)
    return array


@lru_cache(maxsize=10000)
def protein_features(sequence):
    counts = [sequence.count(aa)/len(sequence) for aa in AMINO_ACIDS]
    dipeptides = {a+b:i for i,(a,b) in enumerate(( (a,b) for a in AMINO_ACIDS for b in AMINO_ACIDS))}
    vector = np.zeros(400,dtype=np.float32)
    for index in range(len(sequence)-1):
        code = dipeptides.get(sequence[index:index+2])
        if code is not None:
            vector[code] += 1
    vector /= max(len(sequence)-1,1)
    return np.concatenate([np.asarray(counts,dtype=np.float32),vector,[np.log1p(len(sequence))/10]]).astype(np.float32)


def features(frame):
    return np.asarray([np.concatenate([drug_features(smi),protein_features(seq)]) for smi,seq in zip(frame.smiles,frame.protein_sequence)],dtype=np.float32)


def applicability(train, query):
    """Exact max Morgan Tanimoto and exact 3-mer Jaccard to training entities."""
    drug_pool = [fingerprint(smiles) for smiles in sorted(set(train.smiles))]
    protein_pool = [kmers(sequence) for sequence in sorted(set(train.protein_sequence))]
    drug_values = {smiles:max(DataStructs.BulkTanimotoSimilarity(fingerprint(smiles),drug_pool),default=0)
                   for smiles in set(query.smiles)}
    protein_values = {sequence:max((kmer_jaccard(kmers(sequence),other) for other in protein_pool),default=0)
                      for sequence in set(query.protein_sequence)}
    result = query[[column for column in ("measurement_id","pair_id","compound_id","target_id","smiles","protein_sequence") if column in query]].copy()
    result["max_train_morgan_tanimoto"] = query.smiles.map(drug_values)
    result["max_train_protein_3mer_jaccard"] = query.protein_sequence.map(protein_values)
    result["within_encoder_lengths"] = (query.smiles.str.len()<=256)&(query.protein_sequence.str.len()<=1024)
    # Prespecified operational flag, not a calibrated safety/confidence bound.
    result["applicability_flag"] = (result.max_train_morgan_tanimoto>=0.4)&(result.max_train_protein_3mer_jaccard>=0.3)&result.within_encoder_lengths
    return result
