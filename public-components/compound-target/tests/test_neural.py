import numpy as np
import pytest
torch = pytest.importorskip("torch")
from nutriomics_dti.neural import AffinityCNN,encode


def test_cnn_shapes_gradient_and_deterministic_eval():
    torch.manual_seed(7)
    model = AffinityCNN()
    drug = torch.tensor(np.stack([encode("CCO",32),encode("CCN",32)])).long()
    protein = torch.tensor(np.stack([encode("ACDEFGHIKLMNPQRSTVWY",64,True)]*2)).long()
    model.train()
    output = model(drug,protein)
    assert output.shape==(2,)
    output.square().mean().backward()
    assert model.drug.embedding.weight.grad.abs().sum()>0
    assert model.protein.embedding.weight.grad.abs().sum()>0
    model.eval()
    assert torch.equal(model(drug,protein),model(drug,protein))


def test_padding_and_unknown_encodings():
    assert encode("C",4).tolist()==[ord("C"),0,0,0]
    assert encode("?",4,True).tolist()==[1,0,0,0]
