from itertools import permutations

import pytest

from schubmult import Permutation


def all_perms(n):
    return [Permutation(p) for p in permutations(range(1, n + 1))]


def bruhat_closure(n):
    """Bruhat order as the transitive closure of the length-increasing transposition covers."""
    perms = all_perms(n)
    covers = {p: set() for p in perms}
    for p in perms:
        for i in range(n):
            for j in range(i + 1, n):
                q = p.swap(i, j)
                if q.inv == p.inv + 1:
                    covers[p].add(q)
    reach = {p: {p} for p in perms}
    changed = True
    while changed:
        changed = False
        for p in sorted(perms, key=lambda x: -x.inv):
            for q in covers[p]:
                if not reach[q] <= reach[p]:
                    reach[p] |= reach[q]
                    changed = True
    return reach


# ---------------------------------------------------------------- bruhat_leq


def test_bruhat_leq_basic():
    identity = Permutation([])
    assert identity.bruhat_leq(identity)
    assert identity.bruhat_leq(Permutation([3, 1, 2]))
    assert Permutation([2, 1, 3]).bruhat_leq(Permutation([3, 1, 2]))
    assert not Permutation([3, 1, 2]).bruhat_leq(Permutation([2, 1, 3]))
    assert Permutation([2, 1]).bruhat_leq(Permutation([2, 1, 3]))


def test_bruhat_leq_reflexive_and_antisymmetric():
    for u in all_perms(4):
        assert u.bruhat_leq(u)
        for w in all_perms(4):
            if u.bruhat_leq(w) and w.bruhat_leq(u):
                assert u == w


def test_bruhat_leq_requires_lex_is_not_enough():
    """Regression: sorted-prefix comparison must be componentwise, not lexicographic.

    [1, 4] is lexicographically <= [2, 3] but 4 > 3, so [1,4,2,3] is not below
    [2,3,4,1] in Bruhat order.
    """
    assert not Permutation([1, 4, 2, 3]).bruhat_leq(Permutation([2, 3, 4, 1]))
    assert not Permutation([1, 4, 2, 3]).bruhat_leq(Permutation([3, 2, 4, 1]))
    assert not Permutation([1, 4, 3, 2]).bruhat_leq(Permutation([3, 2, 4, 1]))


def test_bruhat_leq_length_monotone():
    for u in all_perms(4):
        for w in all_perms(4):
            if u.bruhat_leq(w):
                assert u.inv <= w.inv


def test_bruhat_leq_bounds():
    n = 4
    identity = Permutation([])
    w0 = Permutation.w0(n)
    for u in all_perms(n):
        assert identity.bruhat_leq(u)
        assert u.bruhat_leq(w0)


@pytest.mark.parametrize("n", [3, 4, 5])
def test_bruhat_leq_matches_cover_closure(n):
    reach = bruhat_closure(n)
    for u in all_perms(n):
        for w in all_perms(n):
            assert u.bruhat_leq(w) == (w in reach[u]), f"{list(u)} <= {list(w)}"


def test_bruhat_leq_transitive():
    perms = all_perms(4)
    for u in perms:
        for v in perms:
            if not u.bruhat_leq(v):
                continue
            for w in perms:
                if v.bruhat_leq(w):
                    assert u.bruhat_leq(w)


def test_bruhat_leq_inverse_invariant():
    for u in all_perms(4):
        for w in all_perms(4):
            assert u.bruhat_leq(w) == (~u).bruhat_leq(~w)


# ------------------------------------------------------------ weak_order_leq


def test_weak_order_leq_basic():
    identity = Permutation([])
    assert identity.weak_order_leq(identity)
    assert identity.weak_order_leq(Permutation([3, 1, 2]))
    assert Permutation([2, 1, 3]).weak_order_leq(Permutation([3, 1, 2]))
    assert not Permutation([3, 1, 2]).weak_order_leq(Permutation([2, 1, 3]))


def test_weak_order_leq_is_inversion_containment():
    for u in all_perms(4):
        for w in all_perms(4):
            assert u.weak_order_leq(w) == u.inversion_set.issubset(w.inversion_set)


def test_weak_order_leq_refines_bruhat():
    for u in all_perms(4):
        for w in all_perms(4):
            if u.weak_order_leq(w):
                assert u.bruhat_leq(w)


def test_weak_order_leq_is_length_additive():
    """u <= w in weak order iff length(w u^-1) = length(w) - length(u)."""
    for u in all_perms(4):
        for w in all_perms(4):
            assert u.weak_order_leq(w) == ((w * (~u)).inv == w.inv - u.inv)


def test_weak_order_leq_partial_order():
    perms = all_perms(4)
    for u in perms:
        assert u.weak_order_leq(u)
        for w in perms:
            if u.weak_order_leq(w) and w.weak_order_leq(u):
                assert u == w
            if not u.weak_order_leq(w):
                continue
            for v in perms:
                if w.weak_order_leq(v):
                    assert u.weak_order_leq(v)


def test_weak_order_leq_bounds():
    n = 4
    identity = Permutation([])
    w0 = Permutation.w0(n)
    for u in all_perms(n):
        assert identity.weak_order_leq(u)
        assert u.weak_order_leq(w0)


def test_weak_order_leq_is_strictly_weaker_than_bruhat():
    u = Permutation([1, 3, 2])
    w = Permutation([3, 1, 2])
    assert u.bruhat_leq(w)
    assert not u.weak_order_leq(w)


# ----------------------------------------------------------- weak_order_meet


def test_weak_order_meet_basic():
    identity = Permutation([])
    u = Permutation([3, 1, 2])
    assert u.weak_order_meet(u) == u
    assert u.weak_order_meet(identity) == identity
    assert identity.weak_order_meet(u) == identity
    assert Permutation([2, 1, 3]).weak_order_meet(Permutation([3, 1, 2])) == Permutation([2, 1, 3])


def test_weak_order_meet_of_incomparable():
    u = Permutation([3, 1, 2])
    w = Permutation([1, 3, 2])
    assert not u.weak_order_leq(w)
    assert not w.weak_order_leq(u)
    assert u.weak_order_meet(w) == Permutation([])


def test_weak_order_meet_commutative():
    for u in all_perms(4):
        for w in all_perms(4):
            assert u.weak_order_meet(w) == w.weak_order_meet(u)


def test_weak_order_meet_idempotent_and_absorbing():
    for u in all_perms(4):
        assert u.weak_order_meet(u) == u
        assert u.weak_order_meet(Permutation.w0(4)) == u


@pytest.mark.parametrize("n", [3, 4])
def test_weak_order_meet_is_greatest_lower_bound(n):
    perms = all_perms(n)
    for u in perms:
        for w in perms:
            m = u.weak_order_meet(w)
            assert m.weak_order_leq(u)
            assert m.weak_order_leq(w)
            for x in perms:
                if x.weak_order_leq(u) and x.weak_order_leq(w):
                    assert x.weak_order_leq(m)


def test_weak_order_meet_associative():
    perms = all_perms(4)
    for u in perms:
        for v in perms:
            for w in perms:
                left = u.weak_order_meet(v).weak_order_meet(w)
                right = u.weak_order_meet(v.weak_order_meet(w))
                assert left == right


# ----------------------------------------------------------- weak_order_join


def test_weak_order_join_basic():
    identity = Permutation([])
    u = Permutation([3, 1, 2])
    assert u.weak_order_join(u) == u
    assert u.weak_order_join(identity) == u
    assert identity.weak_order_join(u) == u
    assert Permutation([2, 1, 3]).weak_order_join(Permutation([3, 1, 2])) == Permutation([3, 1, 2])


def test_weak_order_join_of_incomparable():
    u = Permutation([3, 1, 2])
    w = Permutation([1, 3, 2])
    assert u.weak_order_join(w) == Permutation.w0(3)


def test_weak_order_join_commutative():
    for u in all_perms(4):
        for w in all_perms(4):
            assert u.weak_order_join(w) == w.weak_order_join(u)


def test_weak_order_join_idempotent_and_absorbing():
    for u in all_perms(4):
        assert u.weak_order_join(u) == u
        assert u.weak_order_join(Permutation([])) == u


@pytest.mark.parametrize("n", [3, 4])
def test_weak_order_join_is_least_upper_bound(n):
    perms = all_perms(n)
    for u in perms:
        for w in perms:
            j = u.weak_order_join(w)
            assert u.weak_order_leq(j)
            assert w.weak_order_leq(j)
            for x in perms:
                if u.weak_order_leq(x) and w.weak_order_leq(x):
                    assert j.weak_order_leq(x)


def test_weak_order_join_associative():
    perms = all_perms(4)
    for u in perms:
        for v in perms:
            for w in perms:
                left = u.weak_order_join(v).weak_order_join(w)
                right = u.weak_order_join(v.weak_order_join(w))
                assert left == right


# ------------------------------------------------------------- meet vs join


def test_weak_order_absorption_laws():
    perms = all_perms(4)
    for u in perms:
        for w in perms:
            assert u.weak_order_meet(u.weak_order_join(w)) == u
            assert u.weak_order_join(u.weak_order_meet(w)) == u


def test_weak_order_meet_join_consistent_with_leq():
    perms = all_perms(4)
    for u in perms:
        for w in perms:
            leq = u.weak_order_leq(w)
            assert (u.weak_order_meet(w) == u) == leq
            assert (u.weak_order_join(w) == w) == leq


def test_weak_order_join_is_w0_dual_of_meet():
    n = 4
    w0 = Permutation.w0(n)
    for u in all_perms(n):
        for w in all_perms(n):
            dual = (u * w0).weak_order_meet(w * w0) * w0
            assert u.weak_order_join(w) == dual
