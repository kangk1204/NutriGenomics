"""Fresh DearDTI-derived GINE / amino-acid CNN / fusion native-pKd head."""
from __future__ import annotations
from functools import lru_cache
import numpy as np
from rdkit import Chem
import torch
from torch import nn
from .vendor.deardti_encoders import DrugGNNEncoder, AACNNEncoder
from .vendor.deardti_fusion import CrossAttentionFusion, BilinearInteraction

ATOM_DIM = 11
BOND_DIM = 10


@lru_cache(maxsize=40000)
def graph(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not 1 <= molecule.GetNumAtoms() <= 128:
        raise ValueError("Graph requires 1..128 atoms; never silently truncate molecular structure")
    nodes = []
    for atom in molecule.GetAtoms():
        nodes.append([atom.GetAtomicNum()/100, atom.GetDegree()/6, atom.GetFormalCharge()/4,
                      atom.GetTotalNumHs()/4, atom.GetIsAromatic(), atom.IsInRing(),
                      float(atom.GetChiralTag())/4, float(atom.GetHybridization())/8,
                      atom.GetMass()/200, atom.GetNumRadicalElectrons()/2, atom.GetIsotope()/250])
    edges, attributes = [], []
    for bond in molecule.GetBonds():
        kind = bond.GetBondType()
        values = [float(kind==type_) for type_ in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)]
        values += [float(bond.GetIsConjugated()), float(bond.IsInRing()), float(bond.GetStereo())/6,
                   float(bond.GetBondDir())/6, bond.GetBondTypeAsDouble()/3, float(kind==Chem.BondType.UNSPECIFIED)]
        for start,end in ((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()),(bond.GetEndAtomIdx(),bond.GetBeginAtomIdx())):
            edges.append((start,end)); attributes.append(values)
    return np.asarray(nodes,np.float32), np.asarray(edges,np.int64).reshape(-1,2), np.asarray(attributes,np.float32).reshape(-1,BOND_DIM)


class PairDataset(torch.utils.data.Dataset):
    def __init__(self, frame, mean=0., scale=1.):
        self.frame = frame.reset_index(drop=True)
        self.mean, self.scale = mean, scale
        self.graphs = {value: graph(value) for value in set(frame.smiles)}
    def __len__(self):
        return len(self.frame)
    def __getitem__(self, index):
        row = self.frame.iloc[index]
        return self.graphs[row.smiles], row.protein_sequence, (float(row.get("pkd",self.mean))-self.mean)/self.scale


def collate(items):
    nmax = max(len(item[0][0]) for item in items)
    batch = len(items)
    nodes = torch.zeros(batch,nmax,ATOM_DIM)
    adj = torch.zeros(batch,nmax,nmax)
    edges = torch.zeros(batch,nmax,nmax,BOND_DIM)
    mask = torch.zeros(batch,nmax)
    labels,seqs = [],[]
    for i,((atoms,index,attrs),sequence,label) in enumerate(items):
        n = len(atoms)
        nodes[i,:n] = torch.from_numpy(atoms)
        mask[i,:n] = 1
        if len(index):
            adj[i,index[:,0],index[:,1]] = 1
            edges[i,index[:,0],index[:,1]] = torch.from_numpy(attrs)
        labels.append(label); seqs.append(sequence)
    return {"nodes":nodes,"adj":adj,"edge_feats":edges,"node_mask":mask,"seqs":seqs,"labels":torch.tensor(labels,dtype=torch.float32)}


class NativeKdGINE(nn.Module):
    def __init__(self, hidden=64, max_tokens=128, protein_cap=1024):
        super().__init__()
        self.drug = DrugGNNEncoder(ATOM_DIM,BOND_DIM,hidden=hidden,n_layers=3,jk=True)
        self.target = AACNNEncoder(hidden=hidden,emb_dim=64,cap=protein_cap,max_tokens=max_tokens)
        self.fusion = CrossAttentionFusion(hidden=hidden,n_heads=4,dropout=.1)
        self.bilinear = BilinearInteraction(hidden=hidden,rank=16)
        self.head = nn.Sequential(nn.Linear(hidden*3,128),nn.LayerNorm(128),nn.ReLU(),nn.Dropout(.1),nn.Linear(128,64),nn.ReLU(),nn.Linear(64,1))
    def forward(self,batch):
        atoms = self.drug(batch["nodes"],batch["adj"],batch["edge_feats"],batch["node_mask"])
        residues,residue_mask = self.target(batch["seqs"])
        drug,target = self.fusion(atoms,residues,batch["node_mask"],residue_mask)
        joint,_ = self.bilinear(atoms,residues,batch["node_mask"],residue_mask)
        return self.head(torch.cat((drug,target,joint),dim=1)).squeeze(-1)


def predict(model,frame,device,mean,scale,batch_size=32):
    model.eval()
    loader = torch.utils.data.DataLoader(PairDataset(frame,mean,scale),batch_size=batch_size,shuffle=False,collate_fn=collate,num_workers=0)
    predictions = []
    with torch.inference_mode():
        for batch in loader:
            batch = {key:value.to(device) if isinstance(value,torch.Tensor) else value for key,value in batch.items()}
            predictions.extend((model(batch)*scale+mean).cpu().numpy().tolist())
    return np.asarray(predictions)
