"""Tests for ``schubmult.mult.separated_descents`` (pipe-puzzle double Grothendieck product)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import beta, y, z
from schubmult.mult.double import schubmult_double
from schubmult.mult.separated_descents import grothmult_double, separated_descents_coeffs
from schubmult.symbolic import S, sympify_sympy

S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
sp = sympify_sympy


def _separated_pairs(perms):
    for u in perms:
        for v in perms:
            des_u = u.descents(zero_indexed=False)
            des_v = v.descents(zero_indexed=False)
            max_u = max(des_u) if des_u else 0
            min_v = min(des_v) if des_v else None
            if min_v is None or max_u <= min_v:
                yield u, v


def _same(d1, d2):
    for w in set(d1) | set(d2):
        if sympy.cancel(sp(d1.get(w, 0)) - sp(d2.get(w, 0))) != 0:
            return False
    return True


def test_grothmult_double_identity_pairs_always_separated():
    idp = Permutation([])
    for v in S3:
        assert _same(grothmult_double({idp: S.One}, v, z, y, beta), {v: 1})


def test_grothmult_double_reduces_to_schubmult_double_at_beta0():
    # G_u(x, y) G_v(x, t) = sum_w c^w(t, y) G_w(x, t); at beta = 0 this is the classical
    # (commuted) double Schubert product S_v(x, t) S_u(x, y) = sum_w c^w(t, y) S_w(x, t).
    for u, v in _separated_pairs(S4):
        sep = grothmult_double({u: S.One}, v, z, y, S.Zero)
        classical = schubmult_double({v: S.One}, u, z, y)
        assert _same(sep, classical), (u, v)


def test_separated_descents_coeffs_matches_grothmult_double():
    for u, v in _separated_pairs(S3):
        single = separated_descents_coeffs(u, v, z, y, beta)
        full = grothmult_double({u: S.One}, v, z, y, beta)
        assert _same(single, full), (u, v)


def test_separated_descents_coeffs_raises_without_separated_descents():
    import pytest

    # [1, 3, 2] has its only descent at 2, [2, 1, 3] has its only descent at 1: not separated.
    with pytest.raises(ValueError):
        separated_descents_coeffs(Permutation([1, 3, 2]), Permutation([2, 1, 3]), z, y, beta)
