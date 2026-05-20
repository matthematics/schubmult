import itertools


def test_forest_rc_ring_principal_products_n3():
    from schubmult import Permutation, Sx
    from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
    from schubmult.rings.polynomial_algebra import PolynomialAlgebra, ForestPolyBasis

    n = 3
    comp_length = n - 1
    perms = Permutation.all_permutations(n)

    forest_ring = ForestRCGraphRing()
    forest_poly = PolynomialAlgebra(ForestPolyBasis(Sx.genset))

    for perm1, perm2 in itertools.product(perms, repeat=2):
        comp1 = tuple(perm1.pad_code(comp_length))
        comp2 = tuple(perm2.pad_code(comp_length))

        expected = forest_poly(*comp1) * forest_poly(*comp2)
        forest1 = forest_ring.forest_poly(comp1)
        forest2 = forest_ring.forest_poly(comp2)

        principal_sum = 0
        for rc1, coeff1 in forest1.items():
            for rc2, coeff2 in forest2.items():
                product = forest_ring(rc1) * forest_ring(rc2)
                for out_rc, coeff in product.items():
                    if out_rc.is_principal:
                        principal_sum += coeff1 * coeff2 * coeff * forest_poly(*out_rc.forest_weight)

        assert expected.almosteq(principal_sum), (
            f"Forest principal mismatch for {perm1.trimcode} * {perm2.trimcode}: "
            f"expected {expected}, got {principal_sum}"
        )
