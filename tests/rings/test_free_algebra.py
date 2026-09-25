import pytest


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

from schubmult.rings.free_algebra import (
    WordBasis,
    SchubertBasis,
    ElementaryBasis,
    FundamentalSlideBasis,
    ForestBasis,
    MonomialSlideBasis,
    FA,
    KeyBasis,
    SchurElementaryBasis,
    GrothendieckBasis,
    GroveBasis,
    GlideBasis,
    LascouxBasis,
)


@pytest.mark.parametrize("basis", [WordBasis, SchubertBasis, ElementaryBasis, FundamentalSlideBasis, ForestBasis, MonomialSlideBasis, KeyBasis, SchurElementaryBasis, GrothendieckBasis, GroveBasis, GlideBasis, LascouxBasis])
def test_word_basis_transitions(basis):
    word_elem = FA(0, 2, 0, 1)

    word_elem2 = word_elem.change_basis(basis).change_basis(WordBasis)
    assert word_elem2 == word_elem


@pytest.mark.parametrize("basis", [WordBasis, SchubertBasis, ElementaryBasis, FundamentalSlideBasis, ForestBasis, MonomialSlideBasis, KeyBasis, SchurElementaryBasis, GrothendieckBasis, GroveBasis, GlideBasis, LascouxBasis])
def test_schubert_basis_transitions(basis):
    from schubmult import ASx, uncode

    schub_elem = ASx(uncode([0, 3, 0, 1]))

    schub_elem2 = schub_elem.change_basis(basis).change_basis(SchubertBasis)
    assert schub_elem2 == schub_elem


@pytest.mark.parametrize("basis", [WordBasis, SchubertBasis, ElementaryBasis, FundamentalSlideBasis, ForestBasis, MonomialSlideBasis, KeyBasis, SchurElementaryBasis, GrothendieckBasis, GroveBasis, GlideBasis, LascouxBasis])
def test_basis_product_matches_word_concatenation(basis):
    left_word = FA(2)
    right_word = FA(2, 1)
    left = left_word.change_basis(basis)
    right = right_word.change_basis(basis)

    assert (left * right).change_basis(WordBasis) == left_word * right_word


@pytest.mark.parametrize("basis", [WordBasis, SchubertBasis, ElementaryBasis, FundamentalSlideBasis, ForestBasis, MonomialSlideBasis, KeyBasis, SchurElementaryBasis, GrothendieckBasis, GroveBasis, GlideBasis, LascouxBasis])
def test_basis_coproduct_matches_word_coproduct(basis):
    from schubmult.rings.free_algebra.free_algebra_basis import FreeAlgebraBasis

    word_elem = FA(1, 0, 2) + 2 * FA(0, 1) - FA(2, 1)
    cprd = word_elem.change_basis(basis).coproduct()

    assert FreeAlgebraBasis.change_tensor_basis(cprd, WordBasis, WordBasis) == word_elem.coproduct()


def test_word_coproduct_splits_letters():
    T = FA @ FA
    assert FA(2).coproduct() == T(((0,), (2,))) + T(((1,), (1,))) + T(((2,), (0,)))
    assert FA(1, 0).coproduct() == T(((0, 0), (1, 0))) + T(((1, 0), (0, 0)))


def test_word_product_is_concatenation():
    assert FA(1, 0, 2) * FA(2, 1) == FA(1, 0, 2, 2, 1)
    assert (FA(1) + FA(0, 2)) * FA(3) == FA(1, 3) + FA(0, 2, 3)


def test_grothendieck_product_matches_schubert_product():
    from schubmult import FreeAlgebra, uncode

    AGx = FreeAlgebra(GrothendieckBasis)
    left = AGx(uncode([1, 0, 1]), 3)
    right = AGx(uncode([0, 2]), 2)
    assert (left * right).change_basis(SchubertBasis) == left.change_basis(SchubertBasis) * right.change_basis(SchubertBasis)


# --- SchurElementaryBasis tests ---


def test_schur_elementary_is_key():
    assert SchurElementaryBasis.is_key(((1, 2), (0, 1, 3)))
    assert SchurElementaryBasis.is_key(([1], [0, 2]))
    assert not SchurElementaryBasis.is_key((1, 2, 3))
    assert not SchurElementaryBasis.is_key(((1,),))


def test_schur_elementary_as_key():
    assert SchurElementaryBasis.as_key(([1, 2], [0, 1, 3])) == ((1, 2), (0, 1, 3))


def test_schur_elementary_zero_monom():
    assert SchurElementaryBasis.zero_monom == ((), ())


def test_schur_elementary_to_schubert_trivial_partition():
    """When lambd[-1] == 0, should delegate to ElementaryBasis."""
    from schubmult.symbolic import S

    result = SchurElementaryBasis.transition_schubert((1, 0), (0, 0, 0))
    assert len(result) > 0
    assert all(isinstance(k, tuple) and len(k) == 2 for k in result)
    assert all(v != S.Zero for v in result.values())


def test_schur_elementary_to_schubert_nontrivial():
    """When lambd has a nonzero last entry, uses the Grassmannian product path."""
    from schubmult.symbolic import S

    result = SchurElementaryBasis.transition_schubert((1,), (0, 2))
    assert len(result) > 0
    assert all(v != S.Zero for v in result.values())


def test_schur_elementary_roundtrip_via_schubert():
    """SE -> Schubert -> SE should be the identity."""
    from schubmult import FreeAlgebra

    SE = FreeAlgebra(SchurElementaryBasis)
    key = ((1,), (0, 2))
    elem = SE(key)
    roundtrip = elem.change_basis(SchubertBasis).change_basis(SchurElementaryBasis)
    assert roundtrip == elem


def test_schur_elementary_roundtrip_via_word():
    """SE -> Word -> SE should be the identity."""
    from schubmult import FreeAlgebra

    SE = FreeAlgebra(SchurElementaryBasis)
    key = ((1, 0), (0, 0, 0))
    elem = SE(key)
    roundtrip = elem.change_basis(WordBasis).change_basis(SchurElementaryBasis)
    assert roundtrip == elem


def test_schubert_to_schur_elementary():
    """Schubert -> SE should produce valid SE keys and round-trip back."""
    from schubmult import ASx, uncode

    perm = uncode([2, 1])
    se_elem = ASx(perm).change_basis(SchurElementaryBasis)
    assert len(se_elem) > 0
    # All keys should be (tuple, tuple) pairs
    for k in se_elem:
        assert SchurElementaryBasis.is_key(k)
    # Round-trip back
    schub_back = se_elem.change_basis(SchubertBasis)
    assert schub_back == ASx(perm)


def test_schubert_schur_elementary_roundtrip_all_s5():
    from schubmult import ASx, Permutation

    for perm in Permutation.all_permutations(5):
        elem = ASx(perm, 4)
        assert elem.change_basis(SchurElementaryBasis).change_basis(SchubertBasis) == elem, perm


def test_schur_elementary_printing_term():
    """Printing term should use SE prefix."""
    key = ((1, 2), (0, 1, 3))
    pt = SchurElementaryBasis.printing_term(key)
    from sympy import sstr

    s = sstr(pt)
    assert "SE" in s
