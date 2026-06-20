from schubmult.combinatorics.wc_graph import WCGraph
from schubmult.symbolic import S, expand, Symbol
from schubmult.symbolic.poly.schub_poly import grothendieck_poly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet
from schubmult.combinatorics.inversions_tableau import InversionsTableau
from schubmult.abc import x
from schubmult import Permutation
import pytest

@pytest.mark.xfail(reason="InversionsTableau is currently broken, but this test should pass once it's fixed")
def test_wc_biject():
    for w in Permutation.all_permutations(4):
        wc_graphs = set([InversionsTableau.from_wc_graph(wc) for wc in WCGraph.all_wc_graphs(w, 3)])
        for wc in wc_graphs:
            assert wc.is_valid, f"Invalid WC graph tableau generated for {w}"
        #assert all(wc.is_valid for wc in wc_graphs), f"Invalid WC graph tableau generated for {w}"
        itabs = set(InversionsTableau.all_set_valued_inversions_tableaux(w))
        assert len(wc_graphs) == len(itabs), f"Failed for {w}: {len(wc_graphs)=}, {len(itabs)=}\n{wc_graphs.difference(itabs)}\n{[tab.perm_word for tab in itabs.difference(wc_graphs)]}"

@pytest.mark.xfail(reason="InversionsTableau is currently broken, but this test should pass once it's fixed")
def test_groth_match():
    beta = Symbol("β")
    zz = ZeroGeneratingSet()

    for w in Permutation.all_permutations(4):
        result = S.Zero
        for it in InversionsTableau.all_set_valued_inversions_tableaux(w):
            result += it.polyvalue(x, beta=beta, prop_beta=True)
        via_it = result
        direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        assert expand(via_it - direct, deep=True) == S.Zero, f"Failed for {w}: {expand(via_it - direct, deep=True)=}"
