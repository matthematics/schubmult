def test_basic_forest():
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra
    from schubmult.abc import x
    
    Forest = PolynomialAlgebra(ForestPolyBasis(x))
    assert (Forest(0,2,0,1).expand() - (x[2]**2*x[4] + x[1]*x[2]*x[4] + x[1]**2*x[4] + x[2]**2*x[3] + x[1]*x[2]*x[3] + x[1]**2*x[3] + x[1]**2*x[2] + x[1]*x[2]**2)).expand() == 0

def test_forest_expansion():
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, SchubertPolyBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra
    n = 5
    perms = Permutation.all_permutations(n)
    Sch = PolynomialAlgebra(SchubertPolyBasis(Sx))
    for perm in perms:
        assert (Sch(perm).change_basis(ForestPolyBasis(Sx.genset)).expand() - Sch(perm).expand()).expand() == 0

