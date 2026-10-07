"""Selected DearDTI MIT components adapted for exact native-Kd regression.
Frozen reference cc741b8adfa445e45dedfcef884567cb125a819a.
Changes: GINE normalization sees valid atoms only; AA-CNN masks after every
convolution and normalizes adaptive pooling by valid residue count.
No old labels, weights, KIBA transforms, retrieval or head are imported.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F
_AA_VOCAB = "ACDEFGHIKLMNPQRSTVWYXBZUO"
_AA2I = {a: i + 1 for i, a in enumerate(_AA_VOCAB)}

class GINELayer(nn.Module):
    def __init__(self, hidden):
        super().__init__()
        self.mlp = nn.Sequential(nn.Linear(hidden, hidden), nn.ReLU(),
                                 nn.Linear(hidden, hidden))
        self.eps = nn.Parameter(torch.zeros(1))
        self.bn = nn.BatchNorm1d(hidden)

    def forward(self, x, adj, edge_h, mask):
        B, N, H = x.shape
        msg = F.relu(x.unsqueeze(1) + edge_h)          # [B,N,N,H]
        msg = (msg * adj.unsqueeze(-1)).sum(dim=2)     # [B,N,H]
        out = self.mlp((1 + self.eps) * x + msg)
        valid = mask.bool()
        normalized = torch.zeros_like(out)
        if self.training and int(valid.sum()) == 1:
            normalized[valid] = F.batch_norm(out[valid], self.bn.running_mean, self.bn.running_var, self.bn.weight, self.bn.bias, training=False, eps=self.bn.eps)
        else:
            normalized[valid] = self.bn(out[valid])
        out = normalized
        return F.relu(out) * mask.unsqueeze(-1)

class DrugGNNEncoder(nn.Module):
    def __init__(self, atom_fdim, bond_fdim, hidden=128, n_layers=3, jk=False):
        super().__init__()
        self.node_proj = nn.Linear(atom_fdim, hidden)
        self.edge_proj = nn.Linear(bond_fdim, hidden)
        self.layers = nn.ModuleList([GINELayer(hidden) for _ in range(n_layers)])
        self.jk = jk
        if jk:
            # jumping-knowledge: fuse every layer's atom representation, so both
            # local (few-hop) and global (many-hop) substructure signal reach the head.
            self.jk_proj = nn.Linear(hidden * (n_layers + 1), hidden)

    def forward(self, nodes, adj, edge_feats, mask):
        h = self.node_proj(nodes) * mask.unsqueeze(-1)
        edge_h = self.edge_proj(edge_feats)
        outs = [h]
        for layer in self.layers:
            h = h + layer(h, adj, edge_h, mask)
            outs.append(h)
        if self.jk:
            h = self.jk_proj(torch.cat(outs, dim=-1)) * mask.unsqueeze(-1)
        return h

class AACNNEncoder(nn.Module):
    """Learned amino-acid embedding + CNN over the raw sequence (consumes AA identity).
    Returns (h[B,L<=cap,hidden], mask[B,L]) to match the PLM encoder interface."""
    def __init__(self, hidden=128, emb_dim=128, cap=1000, max_tokens=200):
        super().__init__()
        self.emb = nn.Embedding(len(_AA_VOCAB) + 1, emb_dim, padding_idx=0)
        self.conv = nn.Sequential(
            nn.Conv1d(emb_dim, hidden, 5, padding=2), nn.ReLU(),
            nn.Conv1d(hidden, hidden, 5, padding=2), nn.ReLU(),
            nn.Conv1d(hidden, hidden, 5, padding=2), nn.ReLU())
        self.cap = cap; self.max_tokens = max_tokens
        self._tokcache = {}
    def _tok(self, s):
        t = self._tokcache.get(s)
        if t is None:
            t = [_AA2I.get(c, _AA2I["X"]) for c in s[:self.cap]]
            self._tokcache[s] = t
        return t
    def forward(self, seqs):
        dev = self.emb.weight.device
        # within-batch dedup: a KIBA batch of 256 pairs has ~100 unique targets;
        # conv the unique sequences once, then scatter back (same trick as LoRA ESM dedup)
        uniq, inv = [], []
        seen = {}
        for s in seqs:
            j = seen.get(s)
            if j is None:
                j = len(uniq); seen[s] = j; uniq.append(s)
            inv.append(j)
        toks = [self._tok(s) for s in uniq]
        L = max(len(t) for t in toks)
        x = torch.zeros(len(toks), L, dtype=torch.long, device=dev)
        for i, t in enumerate(toks): x[i, :len(t)] = torch.tensor(t, device=dev)
        mask = (x != 0).float()
        h = self.emb(x)                                  # [B,L,emb]
        h = h.transpose(1, 2)
        for layer in self.conv:
            h = layer(h) * mask.unsqueeze(1)
        h = h.transpose(1, 2)
        h = h * mask.unsqueeze(-1)
        # Pool each real sequence independently: padding from another longer
        # target must not move a short target's adaptive-pooling boundaries.
        pooled = []
        for i, token in enumerate(toks):
            current = h[i, :len(token)]
            if len(token) > self.max_tokens:
                current = F.adaptive_avg_pool1d(current.transpose(0, 1).unsqueeze(0), self.max_tokens).squeeze(0).transpose(0, 1)
            pooled.append(current)
        width = max(len(current) for current in pooled)
        h = torch.zeros(len(pooled), width, h.shape[-1], device=dev, dtype=h.dtype)
        mask = torch.zeros(len(pooled), width, device=dev, dtype=h.dtype)
        for i, current in enumerate(pooled):
            h[i, :len(current)] = current
            mask[i, :len(current)] = 1
        # scatter unique-target embeddings back to the full batch order
        idx = torch.tensor(inv, dtype=torch.long, device=dev)
        return h.index_select(0, idx), mask.index_select(0, idx)
