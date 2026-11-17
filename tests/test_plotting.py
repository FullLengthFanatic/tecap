from tecap.plotting import _basecomp_legend_handles


def test_basecomp_legend_names_the_bands_only():
    labels = [h.get_label() for h in _basecomp_legend_handles()]
    assert labels == ["30% to 50% A", "60% A or more"]
