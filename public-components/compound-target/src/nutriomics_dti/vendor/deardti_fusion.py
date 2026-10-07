"""Fusion: bidirectional cross-attention (atom<->residue) + bilinear interaction map.

Cross-attention lets each atom attend over residues and vice-versa (TransformerCPI /
HyperAttentionDTI lineage). The bilinear map produces an interpretable [N_atom, N_res]
interaction matrix in the DrugBAN style, usable for visualization.

Design reference: design_proposal.md §3.3; roadmap.md M3.
"""
from __future__ import annotations
import torch
import torch.nn as nn
import torch.nn.functional as F


class CrossAttentionFusion(nn.Module):
    def __init__(self, hidden=128, n_heads=8, dropout=0.1):
        super().__init__()
        self.d2p = nn.MultiheadAttention(hidden, n_heads, batch_first=True, dropout=dropout)
        self.p2d = nn.MultiheadAttention(hidden, n_heads, batch_first=True, dropout=dropout)
        self.ln_d = nn.LayerNorm(hidden)
        self.ln_p = nn.LayerNorm(hidden)

    def forward(self, atom_tok, res_tok, atom_mask, res_mask):
        """atom_tok[B,Na,H], res_tok[B,Nr,H], masks[B,N] (1=valid).
        Returns fused drug vec[B,H], fused target vec[B,H]."""
        akp = atom_mask == 0        # key_padding_mask expects True where padded
        rkp = res_mask == 0
        d, _ = self.d2p(atom_tok, res_tok, res_tok, key_padding_mask=rkp)
        p, _ = self.p2d(res_tok, atom_tok, atom_tok, key_padding_mask=akp)
        d = self.ln_d(atom_tok + d)
        p = self.ln_p(res_tok + p)
        # masked mean pool
        dv = (d * atom_mask.unsqueeze(-1)).sum(1) / atom_mask.sum(1, keepdim=True).clamp(min=1)
        pv = (p * res_mask.unsqueeze(-1)).sum(1) / res_mask.sum(1, keepdim=True).clamp(min=1)
        return dv, pv


class BilinearInteraction(nn.Module):
    """DrugBAN-style low-rank bilinear interaction map (interpretable)."""
    def __init__(self, hidden=128, rank=32):
        super().__init__()
        self.U = nn.Linear(hidden, rank, bias=False)
        self.V = nn.Linear(hidden, rank, bias=False)
        self.out = nn.Linear(rank, hidden)

    def forward(self, atom_tok, res_tok, atom_mask, res_mask):
        A = self.U(atom_tok)                              # [B,Na,r]
        R = self.V(res_tok)                               # [B,Nr,r]
        imap = torch.einsum("bik,bjk->bij", A, R)         # [B,Na,Nr] interaction map
        m = atom_mask.unsqueeze(2) * res_mask.unsqueeze(1)
        imap = imap.masked_fill(m == 0, 0.0)
        # pooled joint representation
        joint = torch.einsum("bij,bik->bk", imap, A) / m.sum((1, 2)).clamp(min=1).unsqueeze(-1)
        return self.out(joint), imap
