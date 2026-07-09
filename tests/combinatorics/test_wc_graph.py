from symengine import Symbol, expand
import pytest

from schubmult import Permutation
from schubmult.abc import x
from schubmult.combinatorics.wc_graph import WCGraph
from schubmult.symbolic import S
from schubmult.symbolic.common_polys import grothendieck_poly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet


def test_groth_match():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.symbolic import S, expand, Symbol
    from schubmult.symbolic.common_polys import grothendieck_poly
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet
    from schubmult.abc import x
    from schubmult import Permutation
    
    beta = Symbol("β")
    zz = ZeroGeneratingSet()

    for w in Permutation.all_permutations(4):
        result = S.Zero
        for wc_graph in WCGraph.all_wc_graphs(w, 3):
            result += wc_graph.polyvalue(x, beta=beta, prop_beta=True)
        via_it = result
        direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        assert expand(via_it - direct, deep=True) == S.Zero, f"Failed for {w}: {expand(via_it - direct, deep=True)=}"

def test_groth_transition():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult.symbolic import S, expand, Symbol
    from schubmult.symbolic.common_polys import grothendieck_poly
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet
    from schubmult.abc import x
    from schubmult import Permutation
    
    beta = Symbol("β")
    zz = ZeroGeneratingSet()
    n = 5
    for w in Permutation.all_permutations(n):
        result = S.Zero
        for wc_graph in WCGraph.all_wc_graphs(w, n - 1):
            if len(wc_graph[-1]) > 0:
                continue
            new_wc_graph = wc_graph.zero_out_last_row()
            rc1 = wc_graph._snap_reduced()
            rc2 = new_wc_graph._snap_reduced()
            result += (beta ** (rc2.perm.inv - rc1.perm.inv))* new_wc_graph.polyvalue(x, beta=beta, prop_beta=True)
        via_it = result
        direct = expand(grothendieck_poly(w, x, zz, beta).subs(x[len(wc_graph)], 0), deep=True)
        assert expand(via_it - direct, deep=True) == S.Zero, f"Failed for {w}: {expand(via_it - direct, deep=True)=}\n{via_it=}\n{direct=}"


def test_set_seq_round_trip():
    from schubmult.combinatorics.wc_graph import WCGraph
    from schubmult import Permutation
    n = 5
    for w in Permutation.all_permutations(n):
        for wc_graph in WCGraph.all_wc_graphs(w, n - 1):
            try:
                wc_graph2 = WCGraph.from_reduced_compatible_set_sequence(*wc_graph.to_reduced_compatible_set_sequence(), length=len(wc_graph))
            except Exception as e:
                raise AssertionError(f"{wc_graph.to_reduced_compatible_set_sequence()} did not round trip from {wc_graph.perm_word}") from e
            assert wc_graph == wc_graph2, f"Round trip failed for {w}: got {wc_graph2}, expected {wc_graph}"
        # via_it = result
        # direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        # assert expand(via_it - direct, deep=True) == S.Zero, f"Failed for {w}: {expand(via_it - direct, deep=True)=}"


def test_pull_out_var_hecke_equation_holds():
    w = Permutation([1, 3, 2])
    for k in range(len(w)):
        sols = WCGraph.pull_out_var_hecke(w, k)
        for row, wpp in sols:
            lhs = Permutation.ref_product(*row) @ wpp.shiftup(1)
            assert lhs == w


@pytest.mark.parametrize("n", [3, 4])
def test_grove_wcs_decompose_grothendieck(n):
    from schubmult import Gx, Sx, uncode
    from schubmult.combinatorics.indexed_forests import grove_polynomial

    beta = Gx._beta
    zz = ZeroGeneratingSet()

    for perm in Permutation.all_permutations(n):
        code = perm.pad_code(n - 1)
        all_wcs = WCGraph.all_wc_graphs(uncode(code), n - 1)
        forest_wcs = {wc for wc in all_wcs if wc.forest_weight == wc.length_vector}

        total = S.Zero
        for forest_wc in forest_wcs:
            grove = WCGraph.grove_wcs(forest_wc.forest_weight, n - 1, forest_wc)
            grove_value = sum(wc.polyvalue(Sx.genset, beta=beta, prop_beta=True) for wc in grove)
            expected = beta ** (len(forest_wc.perm_word) - sum(code)) * grove_polynomial(forest_wc.forest_weight, Sx.genset, beta)
            assert (grove_value - expected).expand() == S.Zero, f"Grove mismatch for {forest_wc.forest_weight}"
            total += grove_value

        direct = grothendieck_poly(uncode(code), Sx.genset, zz, beta)
        assert (total - direct).expand() == S.Zero, f"Grothendieck mismatch for {code}"

