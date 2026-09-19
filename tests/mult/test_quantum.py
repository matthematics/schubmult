"""Tests for ``schubmult.mult.quantum`` (quantum single Schubert polynomial multiplication)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.mult import _accel
from schubmult.mult.quantum import _schubmult_q_fast_python, _vars as _qvars, mult_poly_q, schubmult_q, schubmult_q_fast, single_variable
from schubmult.mult.single import schubmult_py
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars
from schubmult.symbolic.poly.variables import GeneratingSet

x = GeneratingSet("x")
q = _vars.q_var
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy


def _same(d1, d2, subs1=None):
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if sympy.expand(a - b) != 0:
            return False
    return True


def test_schubmult_q_identity_v_empty():
    d = {Permutation([2, 1]): S.One, Permutation([1, 3, 2]): 2}
    assert _same(schubmult_q(d, Permutation([])), d)
    assert _same(schubmult_q(d, Permutation([1, 2])), d)


def test_schubmult_q_known_value():
    # S^q_{21} * S^q_{21} = S_{312} + q_1 * S_id (the quantum correction from the extra
    # Bruhat-chain step through the identity).
    res = schubmult_q({Permutation([2, 1]): S.One}, Permutation([2, 1]))
    assert _same(res, {Permutation([3, 1, 2]): 1, Permutation([]): q[1]})


def test_schubmult_q_fast_matches_schubmult_q():
    for u in S3:
        for v in S3:
            assert _same(schubmult_q_fast({u: S.One}, v), schubmult_q({u: S.One}, v)), (u, v)


def test_schubmult_q_reduces_to_schubmult_py_at_q0():
    q0 = {sp(q[i]): 0 for i in range(1, 8)}
    for u in S3:
        for v in S3:
            assert _same(schubmult_q({u: S.One}, v), schubmult_py({u: S.One}, v), subs1=q0), (u, v)


def test_single_variable_matches_schubmult_q_with_s1():
    for u in S3:
        d = {u: S.One}
        assert _same(single_variable(d, 1), schubmult_q(d, Permutation([2, 1])))


def test_mult_poly_q_matches_single_variable_chain():
    for u in S3:
        d = {u: S.One}
        direct = mult_poly_q(d, x[1] * x[2])
        expected = single_variable(single_variable(d, 2), 1)
        assert _same(direct, expected), u


def test_mult_poly_q_pow_add_scalar_and_plain_list_var_x():
    u = Permutation([2, 1])
    d = {u: S.One}
    assert _same(mult_poly_q(d, x[1] ** 2), single_variable(single_variable(d, 1), 1))
    add_res = mult_poly_q(d, x[1] + 3)
    expected_add = {w: c for w, c in single_variable(d, 1).items()}
    for w, c in d.items():
        expected_add[w] = expected_add.get(w, 0) + 3 * c
    assert _same(add_res, expected_add)
    assert mult_poly_q(d, S(5)) == {w: 5 * c for w, c in d.items()}
    assert _same(mult_poly_q(d, x[1], var_x=[x[0], x[1], x[2]]), single_variable(d, 1))


def test_gvars_n_default():
    assert _qvars.n == 100


def test_schubmult_q_fast_pure_python_matches_dispatcher_and_fallback(monkeypatch):
    for u in S3:
        for v in S3:
            assert _same(_schubmult_q_fast_python({u: S.One}, v), schubmult_q_fast({u: S.One}, v)), (u, v)
    monkeypatch.setattr(_accel, "available", False)
    for u in S3:
        for v in S3:
            assert _same(schubmult_q_fast({u: S.One}, v), _schubmult_q_fast_python({u: S.One}, v)), (u, v)


def test_schubmult_q_fast_pure_python_identity_and_s4_cancellation():
    d = {Permutation([2, 1]): S.One}
    assert _schubmult_q_fast_python(d, Permutation([1, 2])) == d
    # a genuine cancelling v-path (sumval == 0) only shows up for S4-sized permutations.
    S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
    for u in S4:
        for v in S4:
            _schubmult_q_fast_python({u: S.One}, v)
