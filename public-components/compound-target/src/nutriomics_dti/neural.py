"""Freshly initialized SMILES/protein CNN; no pretrained or benchmark weights."""
import numpy as np
import torch
from torch import nn

DRUG_LENGTH = 256
PROTEIN_LENGTH = 1024
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYXBZUOJ"
PROTEIN_CODES = {letter:index+2 for index,letter in enumerate(PROTEIN_ALPHABET)}


def encode(text, length, protein=False):
    result = np.zeros(length,dtype=np.uint8)
    for index,letter in enumerate(text[:length]):
        result[index] = PROTEIN_CODES.get(letter,1) if protein else (ord(letter) if 1<ord(letter)<128 else 1)
    return result


class Encoder(nn.Module):
    def __init__(self, vocabulary, embedding=32, channels=64):
        super().__init__()
        self.embedding = nn.Embedding(vocabulary,embedding,padding_idx=0)
        self.convolutions = nn.Sequential(nn.Conv1d(embedding,channels,5,padding=2),nn.ReLU(),
                                          nn.Conv1d(channels,channels,3,padding=1),nn.ReLU())
    def forward(self,tokens):
        mask = tokens.ne(0).unsqueeze(1)
        values = self.convolutions(self.embedding(tokens).transpose(1,2))
        values = values.masked_fill(~mask,-1e4)
        return values.max(dim=2).values


class AffinityCNN(nn.Module):
    def __init__(self):
        super().__init__()
        self.drug = Encoder(128)
        self.protein = Encoder(len(PROTEIN_CODES)+2)
        self.head = nn.Sequential(nn.Linear(128,128),nn.ReLU(),nn.Dropout(0.15),
                                  nn.Linear(128,32),nn.ReLU(),nn.Linear(32,1))
    def forward(self,drug,protein):
        return self.head(torch.cat([self.drug(drug),self.protein(protein)],dim=1)).squeeze(1)


class PairDataset(torch.utils.data.Dataset):
    def __init__(self,frame,mean=0,scale=1):
        self.drugs = np.asarray([encode(text,DRUG_LENGTH) for text in frame.smiles])
        self.proteins = np.asarray([encode(text,PROTEIN_LENGTH,True) for text in frame.protein_sequence])
        self.labels = ((frame.pkd.to_numpy(dtype=np.float32)-mean)/scale) if "pkd" in frame else np.zeros(len(frame),np.float32)
    def __len__(self):
        return len(self.labels)
    def __getitem__(self,index):
        return torch.from_numpy(self.drugs[index]).long(),torch.from_numpy(self.proteins[index]).long(),self.labels[index]


def predict_neural(model,frame,device,mean,scale,batch_size=256):
    dataset = PairDataset(frame)
    loader = torch.utils.data.DataLoader(dataset,batch_size=batch_size,shuffle=False)
    outputs = []
    model.eval()
    with torch.no_grad():
        for drug,protein,_ in loader:
            outputs.append(model(drug.to(device),protein.to(device)).cpu().numpy()*scale+mean)
    return np.concatenate(outputs)
