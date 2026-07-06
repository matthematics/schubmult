from schubmult.combinatorics.wc_graph import WCGraph
from schubmult.symbolic import S, expand, Symbol
from schubmult.symbolic.common_polys import grothendieck_poly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet
from schubmult.combinatorics.inversions_tableau import InversionsTableau
from schubmult.abc import x
from schubmult import Permutation
import pytest

#@pytest.mark.xfail(reason="InversionsTableau is currently broken, but this test should pass once it's fixed")
@pytest.mark.xfail(reason="InversionsTableau is currently broken")
def test_wc_biject():
    for w in Permutation.all_permutations(4):


        invs = InversionsTableau.all_set_valued_inversions_tableaux(w)
        assert all(iv.is_valid for iv in invs), f"1 Invalid inversions tableau generated for {w}"
        graphs = WCGraph.all_wc_graphs(w, 3)
        it_from_graph = set([InversionsTableau.from_wc_graph(wc) for wc in graphs])
        assert len(it_from_graph) == len(set(invs)), f"Failed 2 for {w}: {set(invs).difference(it_from_graph)}\n{[iv.to_wc_graph() for iv in set(invs).difference(it_from_graph)]}"
        inv_wc_graphs = set([InversionsTableau.to_wc_graph(it, 3) for it in invs])
        assert inv_wc_graphs == graphs, f"Failed 3 for {w}: {inv_wc_graphs.difference(graphs)}\n{[tab.perm_word for tab in inv_wc_graphs.difference(graphs)]}"

@pytest.mark.xfail(reason="InversionsTableau is currently broken")
def test_wc_roundtrip():
    for w in Permutation.all_permutations(4):
        graphs = WCGraph.all_wc_graphs(w, 3)
        inv_wc_graphs = set([InversionsTableau.from_wc_graph(wc) for wc in graphs])
        back_graphs = set()
        for wc in inv_wc_graphs:
            assert wc.is_valid, f"Invalid WC graph tableau generated for {w}"
            back_graphs.add(wc.to_wc_graph(length=3))
        assert back_graphs == graphs, f"Failed for {w}: {back_graphs.difference(graphs)}\n{[tab.perm_word for tab in inv_wc_graphs.difference(graphs)]}"
        #assert all(wc.is_valid for wc in wc_graphs), f"Invalid WC graph tableau generated for {w}"
        # itabs = set(InversionsTableau.all_set_valued_inversions_tableaux(w))
        # assert len(wc_graphs) == len(itabs), f"Failed for {w}: {len(wc_graphs)=}, {len(itabs)=}\n{wc_graphs.difference(itabs)}\n{[tab.perm_word for tab in itabs.difference(wc_graphs)]}"

@pytest.mark.xfail(reason="InversionsTableau is currently broken")
def test_groth_match():
    beta = Symbol("β")
    zz = ZeroGeneratingSet()

    for w in Permutation.all_permutations(4):
        result = S.Zero
        for it in InversionsTableau.all_set_valued_inversions_tableaux(w):
            result += it.polyvalue(x, beta=beta, prop_beta=True)
        via_it = result
        direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        assert expand(via_it - direct, deep=True) == S.Zero, f"Failed for {w}: {expand(via_it - direct, deep=True)=}\n{expand(via_it,deep=True)=}\n{direct=}"
