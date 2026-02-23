def test_key_expansion():
    from schubmult.rings.polynomial_algebra import SchubertPolyBasis, KeyPolyBasis
    from schubmult import Sx, Permutation, PolynomialAlgebra
    from schubmult.symbolic import expand
    n = 5
    perms = Permutation.all_permutations(n)
    Sch = PolynomialAlgebra(SchubertPolyBasis(Sx))
    for perm in perms:
        assert expand(Sch(perm).change_basis(KeyPolyBasis(Sx.genset)).expand() - Sch(perm).expand()) == 0