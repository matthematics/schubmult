"""Tests for ``schubmult.mult.quantum_double`` (quantum double Schubert polynomial multiplication)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import x, y, z
from schubmult.mult.double import schubmult_double
from schubmult.mult.quantum import schubmult_q
from schubmult.mult.quantum_double import factor_out_q, mult_poly_q_double, schubmult_q_double, schubmult_q_double_fast, schubpoly_quantum, single_variable
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars
from schubmult.symbolic.poly.variables import ZeroGeneratingSet

q = _vars.q_var
zero = ZeroGeneratingSet()
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy


def _same(d1, d2, subs1=None, subs2=None):
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if subs2:
            b = b.xreplace(subs2)
        if sympy.expand(a - b) != 0:
            return False
    return True


def test_schubmult_q_double_identity_v_empty():
    d = {Permutation([2, 1]): S.One, Permutation([1, 3, 2]): 2}
    assert _same(schubmult_q_double(d, Permutation([]), y, z), d)
    assert _same(schubmult_q_double(d, Permutation([1, 2]), y, z), d)


def test_schubmult_q_double_fast_matches_reference():
    for u in S3:
        for v in S3:
            assert _same(schubmult_q_double_fast({u: S.One}, v, y, z), schubmult_q_double({u: S.One}, v, y, z)), (u, v)


def test_schubmult_q_double_reduces_to_schubmult_double_at_q0():
    q0 = {sp(q[i]): 0 for i in range(1, 8)}
    for u in S3:
        for v in S3:
            assert _same(schubmult_q_double({u: S.One}, v, y, z), schubmult_double({u: S.One}, v, y, z), subs1=q0), (u, v)


def test_schubmult_q_double_reduces_to_schubmult_q_at_yz0():
    yz0 = {sp(y[i]): 0 for i in range(1, 8)} | {sp(z[i]): 0 for i in range(1, 8)}
    for u in S3:
        for v in S3:
            assert _same(schubmult_q_double({u: S.One}, v, y, z), schubmult_q({u: S.One}, v), subs1=yz0), (u, v)


def test_schubpoly_quantum_polynomial_identity():
    # S^q_u(x; y) S^q_v(x; z) = sum_w c^w_{u,v}(y, z) S^q_w(x; y).
    for u in S3:
        for v in S3:
            prod_dict = schubmult_q_double({u: S.One}, v, y, z)
            lhs = sp(schubpoly_quantum(u, x, y, q)) * sp(schubpoly_quantum(v, x, z, q))
            rhs = sum((sp(c) * sp(schubpoly_quantum(w, x, y, q)) for w, c in prod_dict.items()), sympy.Integer(0))
            assert sympy.expand(lhs - rhs) == 0, (u, v)


def test_single_variable_matches_schubmult_q_double_v_zero_alphabet():
    # S^q_{21}(x, 0) = x_1, so multiplying by it is exactly the equivariant quantum Monk rule.
    for u in S3:
        d = {u: S.One}
        assert _same(single_variable(d, 1, y), schubmult_q_double(d, Permutation([2, 1]), y, zero))


def test_mult_poly_q_double_matches_single_variable_chain():
    for u in S3:
        d = {u: S.One}
        direct = mult_poly_q_double(d, x[1] * x[2], x, y)
        expected = single_variable(single_variable(d, 2, y), 1, y)
        assert _same(direct, expected), u


def test_factor_out_q_roundtrip():
    poly = 3 * q[1] ** 2 * y[1] + 5 * q[2] + 7
    q_dict = factor_out_q(poly)
    total = sum((sp(k) * sp(v) for k, v in q_dict.items()), sympy.Integer(0))
    assert sympy.expand(total - sp(poly)) == 0
