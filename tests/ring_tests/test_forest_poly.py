def test_basic_forest():
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra
    from schubmult.abc import x
    
    Forest = PolynomialAlgebra(ForestPolyBasis(x))
    assert (Forest(0,2,0,1).expand() - (x[2]**2*x[4] + x[1]*x[2]*x[4] + x[1]**2*x[4] + x[2]**2*x[3] + x[1]*x[2]*x[3] + x[1]**2*x[3] + x[1]**2*x[2] + x[1]*x[2]**2)).expand() == 0

def test_forest_slide_expansion():
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, FundamentalSlidePolyBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra
    n = 5
    perms = Permutation.all_permutations(n)
    Forest = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    for perm in perms:
        assert (Forest(*perm.code).change_basis(FundamentalSlidePolyBasis(Sx.genset)).expand() - Forest(*perm.code).expand()).expand() == 0

