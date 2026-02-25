def test_schubert_mul():
    from schubmult import ASx, uncode, Sx, Permutation

    perm1 = uncode([2, 0, 2])
    perm2 = uncode([0, 6, 1])
    result = ASx(perm1) * ASx(perm2)
    p, q = 3, 3
    for (k, n), v in result.items():
        k0 = k * ~(perm2.shiftup(p))
        assert k0.inv == k.inv - perm2.inv
        dct = Sx(k0).coproduct(*list(range(1, p + 1)))
        assert dct[(perm1, Permutation([]))] == v


def test_schubert_coprod():
    from schubmult import ASx, uncode, Sx
    from schubmult.symbolic import S

    perm = uncode([2, 3, 0, 2, 1])
    result = ASx(perm).coproduct()
    for ((k1, n1), (k2, n2)), v in result.items():
        dct = Sx(k1) * Sx(k2)
        assert dct.get(perm, S.Zero) == v


def test_schubert_to_word():
    from schubmult import ASx, uncode, WordBasis, Sx
    from schubmult.symbolic import S, prod

    perm = uncode([3, 1, 0, 2, 1])
    result = ASx(perm).change_basis(WordBasis)
    for k, v in result.items():
        dct = Sx(prod([Sx.genset[i + 1] ** k[i] for i in range(len(k))]))
        assert dct.get(perm, S.Zero) == v


def test_word_to_schubert():
    from schubmult import SchubertBasis, Sx, FA
    from schubmult.symbolic import S

    tup = (3, 1, 0, 2, 1)
    word = FA(*tup)
    result = word.change_basis(SchubertBasis)
    for (k, n), v in result.items():
        dct = Sx(k).to_genset_dict()
        assert dct.get(tup, S.Zero) == v


# def test_schubert_to_sepdesc():
#     from schubmult import ASx, uncode, Sx
#     from schubmult.rings.free_algebra import SeparatedDescentsBasis
#     from schubmult.symbolic import S

#     k = 4
#     perm = uncode([3, 1, 0, 2, 1])
#     result = ASx(perm).change_basis(SeparatedDescentsBasis(k))
#     for (k1, k2, n), v in result.items():
#         assert len(k2) <= k
#         assert all(c >= k - 1 for c in k1.descents())
#         dct = Sx(k1) * Sx(k2)
#         assert dct.get(perm, S.Zero) == v


def test_sepdesc_to_schubert():
    from schubmult import uncode, Sx, FreeAlgebra
    from schubmult.rings.free_algebra import SeparatedDescentsBasis, SchubertBasis, WordBasis
    from schubmult.symbolic import S

    k = 4
    ring = FreeAlgebra(SeparatedDescentsBasis(k))
    perm1 = uncode([0, 1, 3, 2, 1])
    perm2 = uncode([0, 2, 1])
    elem = ring(perm1, perm2, 5)
    result = elem.change_basis(SchubertBasis)
    wbelem = elem.change_basis(WordBasis)
    for (k1, n), v in result.items():
        assert n == 5
        poly = Sx(k1).expand()
        assert wbelem.poly_inner_product(poly, Sx.genset, n) == v

def test_schubert_to_elementary():
    from schubmult import ASx, uncode, Sx
    from schubmult.rings.free_algebra import ElementaryBasis, WordBasis
    from schubmult.symbolic import S
    from schubmult.abc import e

    perm = uncode([3, 1, 3, 0, 1])
    result = ASx(perm).change_basis(ElementaryBasis)
    wbelem = ASx(perm).change_basis(WordBasis)
    for (tup, n), v in result.items():
        assert n == 5
        res = Sx([])
        for i, c in enumerate(tup):
            res *= e(c,i+1,res.ring.genset[1:]) if i <= n else e(c,n,res.ring.genset[1:])
        poly = res.expand()
        assert wbelem.poly_inner_product(poly, Sx.genset, n) == v


def test_elementary_to_schubert():
    from schubmult import Sx, FreeAlgebra
    from schubmult.rings.free_algebra import ElementaryBasis, WordBasis, SchubertBasis

    EE = FreeAlgebra(ElementaryBasis)
    tup = (1, 0, 1, 2, 3)
    numvars = 3
    result = EE(tup,numvars).change_basis(SchubertBasis)
    wbelem = EE(tup,numvars).change_basis(WordBasis)
    for (perm, n), v in result.items():
        assert n == numvars
        res = Sx(perm)
        poly = res.expand()
        assert wbelem.poly_inner_product(poly, Sx.genset, n) == v

def test_basis_change_consistency():
    from schubmult import FreeAlgebra, uncode
    from schubmult.rings.free_algebra import ElementaryBasis, WordBasis, SchubertBasis, SeparatedDescentsBasis, ForestBasis, KeyBasis, FundamentalSlideBasis, MonomialSlideBasis

    # bases_to_test = [ElementaryBasis, SchubertBasis, WordBasis, ForestBasis, KeyBasis, FundamentalSlideBasis]
    bases_to_test = [WordBasis,   SchubertBasis,  ForestBasis, KeyBasis, ElementaryBasis, ]
    test_elements = {ElementaryBasis: ((1, 0, 1, 2), 3), SchubertBasis: (uncode([2,0,1,2]), 3), WordBasis: (1,2,0,3), ForestBasis: (2,0,2,1), KeyBasis: (2,0,3)}

    for basis1 in bases_to_test:
        b1 = FreeAlgebra(basis1)
        for basis2 in bases_to_test:
            if basis1 == basis2:
                continue
            b2 = FreeAlgebra(basis2)
            elem = test_elements[basis1]
            result1 = b1(*elem).change_basis(basis2).change_basis(basis1)
            result2 = b1(*elem)
            assert result1 == result2, f"Basis change inconsistency between {basis1.__name__} and {basis2.__name__} for element {elem}"


