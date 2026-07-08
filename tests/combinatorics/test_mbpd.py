"""Tests for the MBPD <-> RCP <-> WCGraph bijection (writing/mbpd.solve.tex)."""

import itertools

import pytest

from schubmult.combinatorics.mbpd import MBPD, RCP
from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.wc_graph import WCGraph


def _all_rcps(n):
    """Every RCP(n), enumerated via all WCGraphs of every w in S_n."""
    for pl in itertools.permutations(range(1, n + 1)):
        w = Permutation(list(pl))
        for wc in WCGraph.all_wc_graphs(w):
            yield w, wc, RCP.from_wcgraph(wc, n=n)


@pytest.mark.parametrize("n", [3, 4])
def test_rothe_valid_and_perm(n):
    for pl in itertools.permutations(range(1, n + 1)):
        w = Permutation(list(pl))
        D = MBPD.rothe(w, n=n)
        assert D.is_valid(), f"Rothe MBPD invalid for {w}"
        assert D.perm() == w, f"Rothe perm mismatch for {w}: got {list(D.perm())}"


@pytest.mark.parametrize("n", [3, 4])
def test_rcp_wcgraph_roundtrip(n):
    for w, wc, rcp in _all_rcps(n):
        assert rcp.is_valid()
        assert rcp.perm() == w
        back = rcp.to_wcgraph()
        assert back.perm_word == wc.perm_word
        assert back.compatible_sequence == wc.compatible_sequence


@pytest.mark.parametrize("n", [3, 4])
def test_phi_psi_are_inverse(n):
    """Psi = Phi^{-1} as a bijection RCP(n) <-> MBPD(n)."""
    mbpd_seen = {}
    for w, _wc, rcp in _all_rcps(n):
        D = rcp.psi()
        assert D.is_valid(), f"Psi(rcp) invalid for {rcp.biletters}"
        # Phi o Psi = id
        assert D.phi().biletters == rcp.biletters
        # weight and permutation preserved
        assert D.perm() == w
        assert tuple(D.weight()) == tuple(rcp.weight())
        # Psi is injective
        key = D._key()
        assert key not in mbpd_seen or mbpd_seen[key] == rcp.biletters
        mbpd_seen[key] = rcp.biletters


@pytest.mark.parametrize("n", [3, 4])
def test_wcgraph_mbpd_bridge(n):
    for w, wc, _rcp in _all_rcps(n):
        D = wc.to_mbpd(n=n)
        assert D.is_valid()
        assert D.perm() == w
        back = WCGraph.from_mbpd(D)
        assert back.perm_word == wc.perm_word
        assert back.compatible_sequence == wc.compatible_sequence


@pytest.mark.parametrize("n", [3, 4])
def test_groth_poly_from_mbpd(n):
    """sum over MBPD_w of beta^(#heavy - ell) x^wt(D) == G^beta_w."""
    from schubmult.abc import x
    from schubmult.symbolic import S, Symbol, expand
    from schubmult.symbolic.common_polys import grothendieck_poly
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet

    beta = Symbol("β")
    zz = ZeroGeneratingSet()

    for pl in itertools.permutations(range(1, n + 1)):
        w = Permutation(list(pl))
        ell = w.inv
        total = S.Zero
        for wc in WCGraph.all_wc_graphs(w):
            D = wc.to_mbpd(n=n)
            assert D.perm() == w
            monomial = S.One
            for i, m in enumerate(D.weight(), start=1):
                if m:
                    monomial *= x[i] ** m
            total += beta ** (D.num_heavy() - ell) * monomial
        direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        assert expand(total - direct, deep=True) == S.Zero, f"Groth mismatch for {w}"
