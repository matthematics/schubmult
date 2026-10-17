"""Tests for probabilistic zero testing in the double and quantum double kernels (``probabilistic=True``):
the shadow evaluator, agreement with the exact kernels after expansion (C++ and pure-Python paths), the
absence of identically-zero coefficients in the output, iterated products with ``q`` already in the input
coefficients, and the CLI flag."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import y, z
from schubmult.mult import _accel
from schubmult.mult._shadow import PRIME, ShadowEvaluator
from schubmult.mult.double import _schubmult_double_python, schubmult_double
from schubmult.mult.quantum_double import _schubmult_q_double_fast_python, schubmult_q_double_fast
from schubmult.symbolic import S, expand, sympify, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars

q = _vars.q_var
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]


def _nonzero_expanded(d):
    """``{w: expanded coefficient}`` without the identically-zero ones: the ground truth both kernels must match."""
    out = {}
    for w, c in d.items():
        e = expand(c)
        if e != 0:
            out[w] = e
    return out


# ---- ShadowEvaluator ---------------------------------------------------------------------------


def test_shadow_evaluator_matches_direct_evaluation():
    ev = ShadowEvaluator(points=2, seed=7)
    expr = sympify((y[1] - z[2]) ** 3 * (y[3] + 2) - sympify(5) / 3 * y[1] * z[1] ** 2 + 11)
    vals = ev(expr)
    assert len(vals) == 2 and all(0 <= v < PRIME for v in vals)
    for i in range(2):
        point = {s: sympify(ev._symbols[s][i]) for s in expr.free_symbols}  # the values the evaluator drew
        direct = sympy.Rational(str(expand(expr.subs(point))))
        assert vals[i] == int(direct.p) * pow(int(direct.q), -1, PRIME) % PRIME


def test_shadow_evaluator_symbols_keep_their_values():
    ev = ShadowEvaluator(points=2, seed=1)
    assert ev(y[1]) == ev(y[1]) == ev(sympify(y[1]))
    a, b = ev(y[1]), ev(y[2])
    assert a != b
    assert ev(y[1] + y[2]) == tuple((u + v) % PRIME for u, v in zip(a, b))
    assert ev(y[1] * y[2]) == tuple(u * v % PRIME for u, v in zip(a, b))


def test_shadow_evaluator_hidden_zero_and_nonzero():
    ev = ShadowEvaluator(points=2, seed=3)
    hidden = sympify((y[1] - z[2]) * (y[3] - z[1]) - (y[1] * y[3] - y[1] * z[1] - z[2] * y[3] + z[1] * z[2]))
    assert expand(hidden) == 0 and hidden != 0  # structurally nonzero, identically zero
    assert ev(hidden) == (0, 0)
    assert any(ev(sympify((y[1] - z[2]) * (y[3] - z[1]) - y[1] * y[3])))


def test_shadow_evaluator_rejects_non_integer_power():
    ev = ShadowEvaluator()
    try:
        ev(sympify(y[1]) ** sympify(1) / 2 ** sympify(y[2]))
    except TypeError:
        return
    raise AssertionError("expected TypeError for a symbolic exponent")


# ---- schubmult_double --------------------------------------------------------------------------


def test_schubmult_double_probabilistic_matches_exact_cpp():
    for u, v in itertools.product(S4, S4):
        for v3 in (y, z):
            exact = _nonzero_expanded(schubmult_double({u: S.One}, v, y, v3))
            prob = schubmult_double({u: S.One}, v, y, v3, probabilistic=True)
            assert _nonzero_expanded(prob) == exact, (u, v)


def test_schubmult_double_probabilistic_matches_exact_python():
    for u, v in itertools.product(S3, S3):
        for v3 in (y, z):
            exact = _nonzero_expanded(_schubmult_double_python({u: S.One}, v, y, v3))
            prob = _schubmult_double_python({u: S.One}, v, y, v3, probabilistic=True)
            assert _nonzero_expanded(prob) == exact, (u, v)


def test_schubmult_double_probabilistic_drops_every_hidden_zero():
    u, v = Permutation([4, 1, 6, 5, 2, 3]), Permutation([3, 6, 1, 5, 2, 4])
    exact = schubmult_double({u: S.One}, v, y, y)
    prob = schubmult_double({u: S.One}, v, y, y, probabilistic=True)
    assert any(expand(c) == 0 for c in exact.values())  # the exact kernel does emit junk here
    assert all(expand(c) != 0 for c in prob.values())
    assert set(prob) == set(_nonzero_expanded(exact))


def test_schubmult_double_probabilistic_with_symbolic_input_coefficients():
    d = {Permutation([2, 1, 3]): y[1] - z[1], Permutation([1, 3, 2]): sympify(2) * y[2]}
    for v in S3:
        assert _nonzero_expanded(schubmult_double(d, v, y, z, probabilistic=True)) == _nonzero_expanded(schubmult_double(d, v, y, z)), v


# ---- schubmult_q_double_fast -------------------------------------------------------------------


def test_schubmult_q_double_fast_probabilistic_matches_exact_cpp():
    for u, v in itertools.product(S4, S4):
        for v3 in (y, z):
            exact = _nonzero_expanded(schubmult_q_double_fast({u: S.One}, v, y, v3, q))
            prob = schubmult_q_double_fast({u: S.One}, v, y, v3, q, probabilistic=True)
            assert _nonzero_expanded(prob) == exact, (u, v)


def test_schubmult_q_double_fast_probabilistic_matches_exact_python():
    for u, v in itertools.product(S3, S3):
        for v3 in (y, z):
            exact = _nonzero_expanded(_schubmult_q_double_fast_python({u: S.One}, v, y, v3, q))
            prob = _schubmult_q_double_fast_python({u: S.One}, v, y, v3, q, probabilistic=True)
            assert _nonzero_expanded(prob) == exact, (u, v)


def test_schubmult_q_double_fast_probabilistic_drops_every_hidden_zero():
    u, v = Permutation([4, 1, 6, 5, 2, 3]), Permutation([3, 6, 1, 5, 2, 4])
    exact = schubmult_q_double_fast({u: S.One}, v, y, y, q)
    prob = schubmult_q_double_fast({u: S.One}, v, y, y, q, probabilistic=True)
    assert any(expand(c) == 0 for c in exact.values())
    assert all(expand(c) != 0 for c in prob.values())
    assert set(prob) == set(_nonzero_expanded(exact))


def test_schubmult_q_double_fast_probabilistic_iterated_product_with_q_in_coefficients():
    # the second factor sees input coefficients that already contain q (kept under the monomial 1)
    u, v, w = Permutation([2, 3, 1]), Permutation([3, 1, 2]), Permutation([1, 3, 2])
    first = schubmult_q_double_fast({u: S.One}, v, y, z, q, probabilistic=True)
    assert any(sympify_sympy(c).has(sympify_sympy(q[1])) for c in first.values())
    exact = schubmult_q_double_fast(schubmult_q_double_fast({u: S.One}, v, y, z, q), w, y, z, q)
    prob = schubmult_q_double_fast(first, w, y, z, q, probabilistic=True)
    assert _nonzero_expanded(prob) == _nonzero_expanded(exact)


def test_schubmult_q_double_fast_probabilistic_reduces_to_double_at_q0():
    q0 = {sympify_sympy(q[i]): 0 for i in range(1, 8)}
    for u, v in itertools.product(S3, S3):
        quantum = schubmult_q_double_fast({u: S.One}, v, y, z, q, probabilistic=True)
        classical = _nonzero_expanded(schubmult_double({u: S.One}, v, y, z, probabilistic=True))
        at_q0 = {}
        for w, c in quantum.items():
            e = sympy.expand(sympify_sympy(c).xreplace(q0))
            if e != 0:
                at_q0[w] = e
        assert {w: sympy.expand(sympify_sympy(c)) for w, c in classical.items()} == at_q0, (u, v)


# ---- _accel wrappers and CLI -------------------------------------------------------------------


def test_accel_wrappers_accept_probabilistic():
    assert _accel.available
    u, v = Permutation([3, 1, 2]), Permutation([2, 3, 1])
    exact = _nonzero_expanded(_accel.schubmult_double({u: S.One}, v, y, z))
    assert _nonzero_expanded(_accel.schubmult_double({u: S.One}, v, y, z, probabilistic=True)) == exact
    exact_q = _nonzero_expanded(_accel.schubmult_q_double_fast({u: S.One}, v, y, z, q))
    assert _nonzero_expanded(_accel.schubmult_q_double_fast({u: S.One}, v, y, z, q, probabilistic=True)) == exact_q


def test_cli_probabilistic_flag():
    from schubmult._scripts.schubmult_double import main as main_double
    from schubmult._scripts.schubmult_q_double import main as main_q_double

    for main in (main_double, main_q_double):
        argv = ["script", "--mixed-var", "--display-mode", "raw", "3", "1", "5", "2", "4", "-", "2", "5", "1", "4", "3"]
        exact = main(argv)
        prob = main([*argv, "--probabilistic"])
        assert isinstance(exact, dict) and isinstance(prob, dict)
        assert _nonzero_expanded({k: sympify(c) for k, c in prob.items()}) == _nonzero_expanded({k: sympify(c) for k, c in exact.items()})
