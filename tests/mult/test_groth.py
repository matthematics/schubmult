"""Tests for ``schubmult.mult.groth`` (single beta-Grothendieck polynomial multiplication)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import beta, x
from schubmult.mult.groth import chevalley_x_k, grothmult_py, mult_poly_groth, single_variable_groth
from schubmult.mult.single import schubmult_py
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import grothendieck_poly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet

zero = ZeroGeneratingSet()
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy


def _groth(v):
    return sp(grothendieck_poly(v, x, zero, beta))


def _same(d1, d2, subs1=None):
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if sympy.expand(a - b) != 0:
            return False
    return True


def test_grothmult_py_identity_v_empty():
    d = {Permutation([2, 1]): S.One, Permutation([1, 3, 2]): 2}
    assert _same(grothmult_py(d, Permutation([])), d)
    assert _same(grothmult_py(d, Permutation([1, 2])), d)


def test_grothmult_py_reduces_to_schubmult_py_at_beta0():
    b0 = {sp(beta): 0}
    for u in S3:
        for v in S3:
            assert _same(grothmult_py({u: S.One}, v), schubmult_py({u: S.One}, v), subs1=b0), (u, v)


def test_grothmult_py_polynomial_identity():
    for u in S3:
        for v in S3:
            prod_dict = grothmult_py({u: S.One}, v)
            rhs = sum((sp(c) * _groth(w) for w, c in prod_dict.items()), sympy.Integer(0))
            assert sympy.expand(_groth(u) * _groth(v) - rhs) == 0, (u, v)


def test_single_variable_groth_matches_grothmult_py_with_s1():
    for u in S3:
        d = {u: S.One}
        assert _same(single_variable_groth(d, 1, beta), grothmult_py(d, Permutation([2, 1])))


def test_mult_poly_groth_matches_single_variable_chain():
    for u in S3:
        d = {u: S.One}
        direct = mult_poly_groth(d, x[1] * x[2], x, beta)
        expected = single_variable_groth(single_variable_groth(d, 2, beta), 1, beta)
        assert _same(direct, expected), u


def test_chevalley_x_k_self_term_cancels():
    for u in S3:
        coeffs = chevalley_x_k(u, 1, beta)
        assert u not in coeffs
