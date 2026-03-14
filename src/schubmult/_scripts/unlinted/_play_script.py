from schubmult import *

if __name__ == "__main__":
    from schubmult.rings.polynomial_algebra.anti_schubert_poly_basis import AntiSchubertPolyBasis
    from schubmult.rings.polynomial_algebra import *

    ASchub = PolynomialAlgebra(AntiSchubertPolyBasis())

    perm = Permutation(uncode([0,5,2]))

    elem = ASchub(perm, 3).change_basis(ASchub._basis.monomial_basis).change_basis(SchubertPolyBasis()).change_basis(KeyPolyBasis(Sx.genset))
    assert Sx.from_expr(elem.expand()) == Sx.from_expr(ASchub(perm,3).expand())
    #change_basis(SchubertPolyBasis())
    print(elem)
    
    
