from symengine import Symbol, expand

from schubmult import Permutation
from schubmult.abc import x
from schubmult.combinatorics.wc_graph import WCGraph
from schubmult.symbolic import S
from schubmult.symbolic.poly.schub_poly import grothendieck_poly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet


def test_wc_graph_grothendieck_matches_s3():
    beta = Symbol("β")
    zz = ZeroGeneratingSet()

    for w in Permutation.all_permutations(3):
        via_wc = expand(WCGraph.grothendieck_polynomial_via_wc(w, x, beta), deep=True)
        direct = expand(grothendieck_poly(w, x, zz, beta), deep=True)
        assert expand(via_wc - direct, deep=True) == S.Zero


def test_pull_out_var_hecke_equation_holds():
    w = Permutation([1, 3, 2])
    for k in range(len(w)):
        sols = WCGraph.pull_out_var_hecke(w, k)
        for row, wpp in sols:
            lhs = Permutation.ref_product(*row) @ wpp.shiftup(1)
            assert lhs == w
