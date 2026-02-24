import pytest

def test_basic_forest():
    from schubmult.rings.polynomial_algebra import Forest
    from schubmult.abc import x
    
    assert (Forest(0,2,0,1).expand() - (x[2]**2*x[4] + x[1]*x[2]*x[4] + x[1]**2*x[4] + x[2]**2*x[3] + x[1]*x[2]*x[3] + x[1]**2*x[3] + x[1]**2*x[2] + x[1]*x[2]**2)).expand() == 0

def test_forest_slide_expansion():
    from schubmult.rings.polynomial_algebra import Forest, FSlide
    from schubmult import Sx, Permutation, PolynomialAlgebra
    n = 5
    perms = Permutation.all_permutations(n)
    for perm in perms:
        assert (Forest(*perm.code).change_basis(FSlide._basis).expand() - Forest(*perm.code).expand()).expand() == 0

def test_schubert_expansion():
    from schubmult.rings.polynomial_algebra import Forest, Schub
    from schubmult import Sx, Permutation, PolynomialAlgebra
    n = 5
    perms = Permutation.all_permutations(n)
    FBasis = Forest._basis
    for perm in perms:
        forestish = Schub(perm).change_basis(FBasis)
        base_expand = Schub(perm).expand()
        assert (forestish.expand().expand() - base_expand.expand()).expand() == 0


def test_forest_product():
    from schubmult.rings.polynomial_algebra import Forest
    
    prod = Forest(0, 2, 0, 1) * Forest(2, 0, 0, 2)
    expected_result = Forest(2, 2, 0, 3) + Forest(2, 2, 1, 2) + Forest(2, 2, 2, 1) + Forest(2, 3, 0, 2) + Forest(2, 3, 1, 1) + Forest(2, 4, 0, 1) + Forest(3, 1, 0, 3) + Forest(3, 1, 1, 2) + Forest(3, 1, 2, 1) + 2*Forest(3, 2, 0, 2) + 2*Forest(3, 2, 1, 1) + 2*Forest(3, 3, 0, 1) + Forest(4, 0, 0, 3) + Forest(4, 0, 1, 2) + Forest(4, 0, 2, 1) + Forest(4, 1, 0, 2) + Forest(4, 1, 1, 1) + 2*Forest(4, 2, 0, 1) + Forest(5, 0, 0, 2) + Forest(5, 0, 1, 1) + Forest(5, 1, 0, 1) + Forest(6, 0, 0, 1)

    assert prod == expected_result
    assert (prod.expand() - expected_result.expand()).expand() == 0
    assert ((Forest(0, 2, 0, 1).expand() * Forest(2, 0, 0, 2).expand()).expand() - expected_result.expand()).expand() == 0


