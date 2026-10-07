import pytest
torch = pytest.importorskip("torch")
from nutriomics_dti.deardti_kd import graph,collate,NativeKdGINE


def test_graph_preserves_stereo_and_rejects_truncation():
    atoms,index,bonds = graph("N[C@@H](C)C(=O)O")
    assert atoms.shape==(6,11) and bonds.shape[1]==10
    assert atoms[:,6].sum()>0
    with pytest.raises(ValueError):
        graph("C"*129)


def test_fresh_native_kd_head_backward_and_padding_mask():
    torch.set_num_threads(1); torch.manual_seed(3)
    sequence="ACDEFGHIKLMNPQRSTVWY"*2
    batch = collate([(graph("CCO"),sequence,0.),(graph("CCN"),sequence+"AA",1.)])
    model = NativeKdGINE()
    prediction = model(batch)
    assert prediction.shape==(2,) and torch.isfinite(prediction).all()
    loss = torch.nn.functional.mse_loss(prediction,batch["labels"])
    loss.backward()
    assert model.drug.node_proj.weight.grad is not None
    assert model.target.emb.weight.grad is not None
    assert model.fusion.d2p.in_proj_weight.grad is not None


def test_no_old_multitask_or_mixed_affinity_head():
    model = NativeKdGINE()
    assert model.head[-1].out_features==1
    assert not hasattr(model,"interaction_head")


def test_encoder_prediction_is_independent_of_padding_from_other_pairs():
    torch.set_num_threads(1); torch.manual_seed(8)
    model = NativeKdGINE().eval()
    short = "ACDEFGHIKLMNPQRSTVWY"*3
    long = "ACDEFGHIKLMNPQRSTVWY"*20
    one = collate([(graph("CCO"),short,0.)])
    padded = collate([(graph("CCO"),short,0.),(graph("C"*20+"O"),long,0.)])
    with torch.inference_mode():
        assert torch.allclose(model(one)[0],model(padded)[0],atol=1e-6)
