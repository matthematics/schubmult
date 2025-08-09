def test_schubert_mul():
    from schubmult import ASx, uncode, Sx

    perm1 = uncode([2, 0, 2])
    perm2 = uncode([0, 6, 1])
    result = ASx(perm1) * ASx(perm2)
    p, q = 3, 3
    for (k, n), v in result.items():
        dct = Sx(k).coproduct(*list(range(1, p + 1)))
        assert dct[(perm1, perm2)] == v


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
    from schubmult import ASx, uncode, SchubertBasis, Sx, FA
    from schubmult.symbolic import S, prod

    tup = (3, 1, 0, 2, 1)
    word = FA(*tup)
    result = word.change_basis(SchubertBasis)
    for (k, n), v in result.items():
        dct = Sx(k).to_genset_dict()
        assert dct.get(tup, S.Zero) == v


def test_schubert_to_sepdesc():
    from schubmult import ASx, uncode, Sx
    from schubmult.rings import SeparatedDescentsBasis
    from schubmult.symbolic import S

    k = 4
    perm = uncode([3, 1, 0, 2, 1])
    result = ASx(perm).change_basis(SeparatedDescentsBasis(k))
    for (k1, k2, n), v in result.items():
        assert len(k2) <= k
        assert all(c >= k - 1 for c in k1.descents())
        dct = Sx(k1) * Sx(k2)
        assert dct.get(perm, S.Zero) == v


def test_sepdesc_to_schubert():
    from schubmult import uncode, Sx
    from schubmult.rings import SeparatedDescentsBasis, SchubertBasis, FreeAlgebra, ASx, WordBasis
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
