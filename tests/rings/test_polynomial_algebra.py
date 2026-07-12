import pytest


def test_basic_forest():
    from schubmult.rings.polynomial_algebra import ForestPoly
    from schubmult.abc import x

    assert (ForestPoly(0, 2, 0, 1).expand() - (x[2] ** 2 * x[4] + x[1] * x[2] * x[4] + x[1] ** 2 * x[4] + x[2] ** 2 * x[3] + x[1] * x[2] * x[3] + x[1] ** 2 * x[3] + x[1] ** 2 * x[2] + x[1] * x[2] ** 2)).expand() == 0


def test_forest_slide_expansion():
    from schubmult.rings.polynomial_algebra import ForestPoly, FSlide
    from schubmult import Permutation

    n = 5
    perms = Permutation.all_permutations(n)
    for perm in perms:
        assert (ForestPoly(*perm.code).change_basis(FSlide._basis).expand() - ForestPoly(*perm.code).expand()).expand() == 0


def test_schubert_expansion():
    from schubmult.rings.polynomial_algebra import ForestPoly, Schub
    from schubmult import Permutation

    n = 5
    perms = Permutation.all_permutations(n)
    fbasis = ForestPoly._basis
    for perm in perms:
        forestish = Schub(perm).change_basis(fbasis)
        base_expand = Schub(perm).expand()
        assert (forestish.expand().expand() - base_expand.expand()).expand() == 0


def test_forest_product():
    from schubmult.rings.polynomial_algebra import ForestPoly

    prod = ForestPoly(0, 2, 0, 1) * ForestPoly(2, 0, 0, 2)
    expected_result = (
        ForestPoly(2, 2, 0, 3)
        + ForestPoly(2, 2, 1, 2)
        + ForestPoly(2, 2, 2, 1)
        + ForestPoly(2, 3, 0, 2)
        + ForestPoly(2, 3, 1, 1)
        + ForestPoly(2, 4, 0, 1)
        + ForestPoly(3, 1, 0, 3)
        + ForestPoly(3, 1, 1, 2)
        + ForestPoly(3, 1, 2, 1)
        + 2 * ForestPoly(3, 2, 0, 2)
        + 2 * ForestPoly(3, 2, 1, 1)
        + 2 * ForestPoly(3, 3, 0, 1)
        + ForestPoly(4, 0, 0, 3)
        + ForestPoly(4, 0, 1, 2)
        + ForestPoly(4, 0, 2, 1)
        + ForestPoly(4, 1, 0, 2)
        + ForestPoly(4, 1, 1, 1)
        + 2 * ForestPoly(4, 2, 0, 1)
        + ForestPoly(5, 0, 0, 2)
        + ForestPoly(5, 0, 1, 1)
        + ForestPoly(5, 1, 0, 1)
        + ForestPoly(6, 0, 0, 1)
    )

    assert prod == expected_result
    assert (prod.expand() - expected_result.expand()).expand() == 0
    assert ((ForestPoly(0, 2, 0, 1).expand() * ForestPoly(2, 0, 0, 2).expand()).expand() - expected_result.expand()).expand() == 0


def test_fundamental_slide_change_basis():
    from schubmult.rings.polynomial_algebra import FundamentalSlidePolyBasis, SchubertPolyBasis, MonomialBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra

    n = 5
    perms = Permutation.all_permutations(n)
    fslide = PolynomialAlgebra(FundamentalSlidePolyBasis(Sx.genset))
    for perm in perms:
        assert (fslide(*perm.code).change_basis(SchubertPolyBasis(Sx)).expand() - fslide(*perm.code).expand()).expand() == 0
        assert (fslide(*perm.code).change_basis(SchubertPolyBasis(Sx)).change_basis(FundamentalSlidePolyBasis(Sx.genset)).expand() - fslide(*perm.code).expand()).expand() == 0
        assert (fslide(*perm.code).change_basis(MonomialBasis(Sx.genset)).expand() - fslide(*perm.code).expand()).expand() == 0


def test_fundamental_slide_product():
    from schubmult.rings.polynomial_algebra import FundamentalSlidePolyBasis, MonomialBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra

    n = 4
    perms = Permutation.all_permutations(n)
    fslide = PolynomialAlgebra(FundamentalSlidePolyBasis(Sx.genset))
    for perm1 in perms:
        for perm2 in perms:
            lhs = (fslide(*perm1.code) * fslide(*perm2.code)).expand().expand()
            rhs = (fslide(*perm1.code).change_basis(MonomialBasis(Sx.genset)).expand() * fslide(*perm2.code).change_basis(MonomialBasis(Sx.genset)).expand()).expand().expand()
            assert (lhs - rhs).expand() == 0


def test_monomial_from_expr():
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra, MonomialBasis
    from schubmult.symbolic import expand
    from schubmult.abc import x

    monom = PolynomialAlgebra(MonomialBasis(x))
    test_poly = (x[1] ** 2 + 2 * x[1] * x[2] + x[2] ** 2) ** 2
    assert expand(monom.from_expr(test_poly).expand() - expand(test_poly)) == 0


def test_monomial_product():
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra, MonomialBasis
    from schubmult.abc import x

    monom = PolynomialAlgebra(MonomialBasis(x))
    poly1 = x[1] ** 2 + 2 * x[1] * x[2] + x[2] ** 2
    poly2 = x[1] + x[2]
    result = poly1 * poly2
    monom1 = monom.from_expr(poly1)
    monom2 = monom.from_expr(poly2)
    result_monom = monom.from_expr(result)
    assert (monom1 * monom2).almosteq(result_monom)
    assert ((monom1 * monom2).expand() - result_monom.expand()).expand() == 0


def test_schubert_key_expansion():
    from schubmult.rings.polynomial_algebra import SchubertPolyBasis, KeyPolyBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra
    from schubmult.symbolic import expand

    n = 5
    perms = Permutation.all_permutations(n)
    sch = PolynomialAlgebra(SchubertPolyBasis(Sx))
    for perm in perms:
        assert expand(sch(perm).change_basis(KeyPolyBasis(Sx.genset)).expand() - sch(perm).expand()) == 0

from schubmult.rings.polynomial_algebra import PolynomialAlgebra, MonomialBasis, SchubertPolyBasis, ElemSymPolyBasis, FundamentalSlidePolyBasis, ForestPolyBasis, MonomialSlidePolyBasis, KeyPolyBasis, GrothendieckPolyBasis, GrovePolyBasis

@pytest.mark.parametrize("basis", [
    MonomialBasis, 
    SchubertPolyBasis, 
    ElemSymPolyBasis, 
    FundamentalSlidePolyBasis, 
    ForestPolyBasis, 
    MonomialSlidePolyBasis, 
    KeyPolyBasis,
    GrothendieckPolyBasis,
    GrovePolyBasis,
    ])
def test_monomial_basis_transitions(basis):
    from schubmult.abc import x

    PA = PolynomialAlgebra(MonomialBasis(x))

    monom_elem = PA(1, 3, 1, 2) + PA(2,2,0,3) + PA(0, 1, 4, 0)
    monom_elem2 = monom_elem.change_basis(basis).change_basis(MonomialBasis(x))
    assert monom_elem2.almosteq(monom_elem)

@pytest.mark.parametrize("basis", [
    MonomialBasis, 
    SchubertPolyBasis, 
    ElemSymPolyBasis, 
    FundamentalSlidePolyBasis, 
    ForestPolyBasis, 
    MonomialSlidePolyBasis, 
    KeyPolyBasis,
    GrothendieckPolyBasis,
    GrovePolyBasis,
    ])
def test_schubert_basis_transitions(basis):
    from schubmult.abc import x
    from schubmult import uncode

    Schub = PolynomialAlgebra(SchubertPolyBasis(x))

    monom_elem = Schub(uncode([2,0,1,3])) - Schub(uncode([2,0,1]), 4)
    monom_elem2 = monom_elem.change_basis(basis).change_basis(SchubertPolyBasis)
    assert monom_elem2.almosteq(monom_elem)
