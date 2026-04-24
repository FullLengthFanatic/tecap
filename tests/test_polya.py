"""Tests for strand-aware polyA atlas handling."""

from tecap.polya import build_polya_index, hits_polya, utr_bin


def test_strand_aware_lookup(synthetic_fixture):
    idx = build_polya_index(synthetic_fixture["polya"], slop=25,
                            allowed_types={"TE", "DS", "DI", "EX"})

    # 9950 on + is in [9900,10000] cluster with PAS hexamer
    hit, has_pas = hits_polya("chr_test", 9950, "+", idx)
    assert hit and has_pas

    # Same position on - is NOT a hit (no cluster there)
    hit, _ = hits_polya("chr_test", 9950, "-", idx)
    assert not hit

    # PAS- cluster: hit but no hexamer
    hit, has_pas = hits_polya("chr_test", 8450, "+", idx)
    assert hit and not has_pas

    # Minus-strand terminal cluster
    hit, has_pas = hits_polya("chr_test", 12050, "-", idx)
    assert hit and has_pas


def test_slop_boundary(synthetic_fixture):
    idx = build_polya_index(synthetic_fixture["polya"], slop=25,
                            allowed_types={"TE", "DS", "DI", "EX"})
    # Cluster [9900,10000] with slop=25 → [9875,10025].
    assert hits_polya("chr_test", 9880, "+", idx)[0]
    assert not hits_polya("chr_test", 9870, "+", idx)[0]
    assert hits_polya("chr_test", 10020, "+", idx)[0]
    assert not hits_polya("chr_test", 10030, "+", idx)[0]


def test_utr_bin_edges():
    assert utr_bin(0)    == 0   # "0-500"
    assert utr_bin(499)  == 0
    assert utr_bin(500)  == 1   # "500-1k"
    assert utr_bin(999)  == 1
    assert utr_bin(1000) == 2   # "1k-2k"
    assert utr_bin(2000) == 3   # "2k-5k"
    assert utr_bin(5000) == 4   # ">5k"
    assert utr_bin(10**9) == 4
