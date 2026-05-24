"""Express forest polynoials in the Schubert basis"""
from schubmult import *
from schubmult.rings.polynomial_algebra import PolynomialAlgebra, ForestPolyBasis, Schub, SchubertPolyBasis, MonomialBasis, PA
from schubmult.symbolic import prod
from schubmult.abc import e, x
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.utils.perm_utils import artin_sequences

def _trim_monomial(index, monom, coeff=1):
    comp1 = [m if i != index + 1 and i != index else 0 if i != index else m - 1 for i, m in enumerate(monom)]
    comp2 = [m if i != index + 1 and i != index else m[index] - 1 if i != index else 0 for i, m in enumerate(monom)]
    dct = {tuple(comp1): coeff}
    dct[tuple(comp2)] = dct.get(tuple(comp2), 0) - coeff
    return dct

def forest_to_schub(comp):
    ForestPoly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    n = max([i + m for i, m in enumerate(comp) if m > 0], default=0)
    if n == 0:
        return Sx([])
    perm = uncode(comp) * Permutation.w0(n + 1)
    monoms = ForestPoly(*perm.pad_code(n)).change_basis(MonomialBasis)
    result = Sx.zero
    for monom in artin_sequences(n):
        weight_vector = [n - i - m for i, m in enumerate(monom)]
        producto = Sx.from_expr(prod([e(w, n - i,x[1:]) for i, w in enumerate(weight_vector) if w > 0]))
        result += PA(*monom).change_basis(ForestPolyBasis).get(pad_tuple(perm.pad_code(n), n), 0) * producto
    return result


if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    for comp in itertools.product(range(n), repeat=n):
        print(f"Forest comp {comp} -> Schub: {forest_to_schub(comp)}")