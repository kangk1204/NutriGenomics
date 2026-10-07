import pytest
from nutriomics_dti.models import metrics


def test_metrics_are_not_mse_and_handle_constant_rank():
    score = metrics([1,2],[2,4])
    assert score["rmse_pkd"]==pytest.approx((2.5)**0.5)
    assert score["mae_pkd"]==pytest.approx(1.5)
    assert score["spearman"]==pytest.approx(1)
    assert metrics([1,2],[3,3])["spearman"] is None
