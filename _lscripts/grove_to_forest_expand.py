from schubmult import *  # noqa
from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor
from schubmult.rings.polynomial_algebra import PolynomialAlgebra
from schubmult.rings.polynomial_algebra.forest_poly_basis import (
    ForestPolyBasis,
    _grove_polynomial_from_indfor,
)
from schubmult.symbolic import Symbol

beta = Gx._beta

ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))


def grove_in_forests(comp):
    forest = weak_composition_to_indfor(comp)
    grove = _grove_polynomial_from_indfor(forest, Sx.genset, beta=beta)
    return ForestPoly.from_expr(grove)


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    perms = Permutation.all_permutations(n)
    comps = sorted({perm.pad_code(n - 1) for perm in perms})
    for comp in comps:
        expansion = grove_in_forests(comp)
        forests = sorted(key for key, coeff in expansion.items() if coeff != 0)
        print(f"grove{comp} -> {forests}")
