def test_change_basis():
    from schubmult.rings.polynomial_algebra import FundamentalSlidePolyBasis, SchubertPolyBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra
    n = 5
    perms = Permutation.all_permutations(n)
    FSlide = PolynomialAlgebra(FundamentalSlidePolyBasis(Sx.genset))
    for perm in perms:
        assert (FSlide(*perm.code).change_basis(SchubertPolyBasis(Sx)).expand() - FSlide(*perm.code).expand()).expand() == 0
        assert (FSlide(*perm.code).change_basis(SchubertPolyBasis(Sx)).change_basis(FundamentalSlidePolyBasis(Sx.genset)).expand() - FSlide(*perm.code).expand()).expand() == 0

def test_product():
    from schubmult.rings.polynomial_algebra import FundamentalSlidePolyBasis, MonomialBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra
    n = 5
    perms = Permutation.all_permutations(n)
    FSlide = PolynomialAlgebra(FundamentalSlidePolyBasis(Sx.genset))
    for perm1 in perms:
        for perm2 in perms:
            assert ((FSlide(*perm1.code)*FSlide(*perm2.code)).expand() - FSlide(*perm1.code).change_basis(MonomialBasis(Sx.genset)) * FSlide(*perm2.code).change_basis(MonomialBasis(Sx.genset)).expand()).expand() == 0