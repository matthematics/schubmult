"""Tests for ``schubmult.mult.transition`` (transition + Monk kernel and the hybrid dispatcher)."""

import itertools
import math
import random

import pytest

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph
from schubmult.mult import _accel
from schubmult.mult.single import schubmult_py
from schubmult.mult.transition import (
    _prefers_transition,
    mult_monomials_py,
    pipe_dream_count,
    schubmult_double_transition,
    schubmult_py_hybrid,
    schubmult_py_transition,
    transition_monomials,
)
from schubmult.symbolic import expand, sympify_sympy
from schubmult.symbolic.poly.schub_poly import schubpoly
from schubmult.symbolic.poly.variables import GeneratingSet, ZeroGeneratingSet

x = GeneratingSet("x")
zero = ZeroGeneratingSet()
S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
S5 = [Permutation(list(p)) for p in itertools.permutations(range(1, 6))]


def _clean(d):
    return {k: v for k, v in d.items() if v != 0}


def _poly_of_monomials(mons):
    return sum(c * math.prod(x[i + 1] ** a for i, a in enumerate(m)) for m, c in mons.items())


@pytest.mark.parametrize("w", S4 + [Permutation([3, 1, 5, 2, 4]), Permutation([4, 5, 1, 3, 2]), Permutation([6, 2, 5, 1, 4, 3])])
def test_transition_monomials_is_schubert_polynomial(w):
    assert expand(sympify_sympy(_poly_of_monomials(transition_monomials(w))) - sympify_sympy(schubpoly(w, x, zero))) == 0


@pytest.mark.parametrize("w", S5)
def test_pipe_dream_count_matches_rc_graphs_and_monomials(w):
    assert pipe_dream_count(w) == RCGraph.count_rc_graphs(w) == sum(transition_monomials(w).values())


@pytest.mark.skipif(not _accel.available, reason="compiled extension")
@pytest.mark.parametrize("w", [Permutation([3, 1, 5, 2, 4]), Permutation([5, 4, 1, 3, 2, 6, 8, 7]), Permutation([2, 7, 1, 5, 3, 8, 4, 6])])
def test_compiled_pipe_dream_count_matches_python(w):
    assert _accel.pipe_dream_count(w) == RCGraph.count_rc_graphs(w)


@pytest.mark.parametrize("u,v", list(itertools.product(S4, S4))[::7])
def test_transition_kernel_matches_vpath_kernel(u, v):
    assert _clean(schubmult_py_transition({u: 1}, v)) == _clean(schubmult_py({u: 1}, v))


def test_python_reference_matches_compiled():
    rng = random.Random(3)
    for _ in range(25):
        n = rng.randint(3, 6)
        u, v = (Permutation(rng.sample(range(1, n + 1), n)) for _ in range(2))
        expected = _clean(schubmult_py({u: 1}, v))
        assert _clean(mult_monomials_py({u: 1}, transition_monomials(v))) == expected
        assert _clean(schubmult_py_transition({u: 1}, v)) == expected
        assert _clean(schubmult_py_hybrid({u: 1}, v)) == expected


def test_linear_combination_and_zero_coefficients():
    d = {Permutation([2, 1, 3]): 2, Permutation([1, 3, 2]): -1, Permutation([3, 1, 2]): 0}
    v = Permutation([2, 3, 1])
    assert _clean(schubmult_py_transition(d, v)) == _clean(schubmult_py(d, v))
    assert _clean(schubmult_py_hybrid(d, v)) == _clean(schubmult_py(d, v))
    assert schubmult_py_hybrid({Permutation([3, 1, 2]): 0}, v) == {}
    assert _clean(schubmult_py_hybrid(d, [])) == _clean(d)


def test_hybrid_matches_on_larger_mixed_products():
    rng = random.Random(11)
    for _ in range(10):
        u, v = (Permutation(rng.sample(range(1, 9), 8)) for _ in range(2))
        assert _clean(schubmult_py_hybrid({u: 1}, v)) == _clean(schubmult_py({u: 1}, v))
    d = {Permutation(rng.sample(range(1, 7), 6)): rng.randint(-3, 3) for _ in range(4)}
    v = Permutation(rng.sample(range(1, 7), 6))
    assert _clean(schubmult_py_hybrid(d, v)) == _clean(schubmult_py(d, v))


def test_cost_model_regimes():
    # a Grassmannian permutation: few v-path transitions, many pipe dreams -> v-path kernel
    assert not _prefers_transition(pd=32768, nvp=69)
    # the "desc2" family: many v-path transitions relative to pipe dreams -> transition kernel
    assert _prefers_transition(pd=2548, nvp=550)
    assert not _prefers_transition(pd=1, nvp=0)


# ---------------------------------------------------------------------------
# double transition
# ---------------------------------------------------------------------------

y = GeneratingSet("y")
z = GeneratingSet("z")


def _same(d1, d2):
    import sympy

    keys = set(d1) | set(d2)
    return all(sympy.expand(sympify_sympy(d1.get(k, 0)) - sympify_sympy(d2.get(k, 0))) == 0 for k in keys)


@pytest.mark.parametrize("w", [Permutation([3, 1, 2]), Permutation([2, 3, 1]), Permutation([3, 4, 1, 2]), Permutation([4, 2, 1, 3]), Permutation([2, 4, 1, 5, 3])])
def test_double_transition_formula(w):
    # S_w(x, z) = (x_r - z_{v(r)}) S_v(x, z) + sum_{i < r, v <. v t_{ir}} S_{v t_{ir}}(x, z)
    import sympy

    from schubmult.mult.transition import _last_descent_split

    r, v = _last_descent_split(w)
    rhs = (x[r] - z[v[r - 1]]) * schubpoly(v, x, z)
    vr, last = v[r - 1], 0
    for i in range(r - 1, 0, -1):
        if last < v[i - 1] < vr:
            last = v[i - 1]
            rhs += schubpoly(v.swap(i - 1, r - 1), x, z)
    assert sympy.expand(sympify_sympy(schubpoly(w, x, z)) - sympify_sympy(rhs)) == 0


@pytest.mark.parametrize("u,v", list(itertools.product(S4, S4))[::11])
def test_double_transition_matches_vpath_kernel(u, v):
    from schubmult.mult.double import schubmult_double

    assert _same(schubmult_double_transition({u: 1}, v, y, z), schubmult_double({u: 1}, v, y, z))


def test_double_transition_python_reference_and_same_alphabet():
    from schubmult.mult.double import schubmult_double

    rng = random.Random(5)
    for _ in range(6):
        n = rng.randint(3, 5)
        u, v = (Permutation(rng.sample(range(1, n + 1), n)) for _ in range(2))
        expected = schubmult_double({u: 1}, v, y, z)
        assert _same(schubmult_double_transition({u: 1}, v, y, z), expected)
        if _accel.available:
            with_python = _accel.schubmult_double_transition
            try:
                _accel.schubmult_double_transition = lambda *a: None
                assert _same(schubmult_double_transition({u: 1}, v, y, z), expected)
            finally:
                _accel.schubmult_double_transition = with_python
        # same alphabet on both sides: the y_{u(r)} - y_{v(r)} scalars vanish on the diagonal
        assert _same(schubmult_double_transition({u: 1}, v, y, y), schubmult_double({u: 1}, v, y, y))


def test_double_transition_linear_combination():
    from schubmult.mult.double import schubmult_double

    d = {Permutation([2, 1, 3]): y[1] - z[2], Permutation([1, 3, 2]): 2, Permutation([3, 1, 2]): 0}
    v = Permutation([2, 3, 1])
    assert _same(schubmult_double_transition(d, v, y, z), schubmult_double(d, v, y, z))
    assert schubmult_double_transition(d, [], y, z) == {k: c for k, c in d.items() if c != 0}
