"""Tests for ``schubmult.mult.groth_double`` (double beta-Grothendieck Monk/Pieri formula)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y, z
from schubmult.mult.double import schubmult_double
from schubmult.mult.groth_double import (
    epsilon_chain,
    grothmult_double,
    grothmult_double_pieri,
    grothmult_double_top,
    monk_chain,
    mult_poly_groth_double,
    single_variable_groth,
)
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import grothendieck_poly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet

zero = ZeroGeneratingSet()
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy
FLIP_YZ = {sp(y[i]): -sp(y[i]) for i in range(1, 12)} | {sp(z[i]): -sp(z[i]) for i in range(1, 12)}


def _groth(v, var2, var3):
    return sp(grothendieck_poly(v, var2, var3, beta))


def _same(d1, d2, subs1=None, subs2=None):
    # coefficients are rational functions in beta/y/z, so compare via cancel, not expand
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if subs2:
            b = b.xreplace(subs2)
        if sympy.cancel(a - b) != 0:
            return False
    return True


def test_grothmult_double_identity_v_empty():
    d = {Permutation([2, 1]): S.One, Permutation([1, 3, 2]): 2}
    assert _same(grothmult_double(d, Permutation([]), y, z, beta), d)
    assert _same(grothmult_double(d, Permutation([1, 2]), y, z, beta), d)


def test_grothmult_double_reduces_to_schubmult_double_at_beta0():
    # G_w(x; y)|_{beta=0} = S_w(x; -y) (see schubmult.mult.groth docstring).
    b0 = {sp(beta): 0}
    for u in S3:
        for v in S3:
            assert _same(grothmult_double({u: S.One}, v, y, z, beta), schubmult_double({u: S.One}, v, y, z), subs1=b0, subs2=FLIP_YZ), (u, v)


def test_grothmult_double_polynomial_identity():
    # G_u(x; y) G_v(x; z) = sum_w c^w_{u,v}(y, z) G_w(x; y), as rational functions.
    for u in S3:
        for v in S3:
            prod_dict = grothmult_double({u: S.One}, v, y, z, beta)
            rhs = sum((sp(c) * _groth(w, x, y) for w, c in prod_dict.items()), sympy.Integer(0))
            assert sympy.cancel(_groth(u, x, y) * _groth(v, x, z) - rhs) == 0, (u, v)


def test_single_variable_groth_double_v_zero_alphabet_matches_top():
    # G_{21}(x, 0) = x_1, so multiplying by it is exactly the equivariant Chevalley rule.
    for u in S3:
        d = {u: S.One}
        assert _same(single_variable_groth(d, 1, y, beta), grothmult_double(d, Permutation([2, 1]), y, zero, beta))


def test_grothmult_double_top_matches_pieri_top_degree():
    # The p = k Pieri term equals the top linear block for k = 1, 2.
    for u in S3:
        for k in (1, 2):
            assert _same(grothmult_double_top({u: S.One}, k, S.Zero, y, beta), grothmult_double_pieri({u: S.One}, k, k, S.Zero, None, y, beta)), (u, k)


def test_mult_poly_groth_double_matches_single_variable_chain():
    for u in S3:
        d = {u: S.One}
        direct = mult_poly_groth_double(d, x[1] * x[2], x, y, beta)
        expected = single_variable_groth(single_variable_groth(d, 2, y, beta), 1, y, beta)
        assert _same(direct, expected), u


def test_monk_chain_is_epsilon_chain_projection():
    for k in (1, 2, 3):
        assert monk_chain(k) == tuple((a, b) for a, b, _ in epsilon_chain(tuple(range(1, k + 1))))


def test_epsilon_chain_single_index_known_shape():
    # A single index k gives (i, k, 0) for i < k and (k, j, 1) for j > k (see docstring).
    k = 2
    chain = epsilon_chain(k, ambient_rank=4)
    for a, b, m in chain:
        if b == k:
            assert m == 0
        elif a == k:
            assert m == 1
