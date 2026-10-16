"""Tests for ``schubmult.utils.perm_utils`` (block lift of arXiv:2610.21053)."""

import itertools

import pytest
import sympy

from schubmult import Permutation
from schubmult.abc import x
from schubmult.mult.single import schubmult_py
from schubmult.symbolic import sympify_sympy
from schubmult.symbolic.poly.schub_poly import schubpoly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet
from schubmult.utils.perm_utils import block_lift, increasing_completion

# Running example of the paper: u = 435126, v = 142635, w = 164325.
U, V, W = Permutation([4, 3, 5, 1, 2, 6]), Permutation([1, 4, 2, 6, 3, 5]), Permutation([1, 6, 4, 3, 2, 5])


def test_increasing_completion():
    assert increasing_completion([1, 4, 5, 6, 8, 9], 9) == Permutation([1, 4, 5, 6, 8, 9, 2, 3, 7])
    assert increasing_completion([3, 1, 2], 5) == Permutation([3, 1, 2])
    assert increasing_completion([], 4) == Permutation([])
    with pytest.raises(ValueError):
        increasing_completion([1, 1], 3)
    with pytest.raises(ValueError):
        increasing_completion([4], 3)


def test_block_lift_example_4_8_and_table_1():
    """Example 4.8 and Table 1 of arXiv:2610.21053 (one-line notation, bar = appended increasing tail)."""
    table = {
        (0, 0, 0, 0, 0, 0): ([1, 4, 2, 6, 3, 5], [1, 6, 4, 3, 2, 5]),
        (1, 0, 0, 1, 0, 0): ([1, 5, 3, 8, 4, 7, 2, 6], [2, 8, 6, 4, 3, 7, 1, 5]),
        (2, 0, 0, 1, 0, 0): ([1, 6, 4, 9, 5, 8, 2, 3, 7], [3, 9, 7, 5, 4, 8, 1, 2, 6]),
        (2, 0, 0, 0, 0, 1): ([1, 6, 4, 8, 5, 7, 2, 3, 9], [3, 9, 6, 5, 4, 7, 1, 2, 8]),
        (3, 0, 0, 1, 0, 0): ([1, 7, 5, 10, 6, 9, 2, 3, 4, 8], [4, 10, 8, 6, 5, 9, 1, 2, 3, 7]),
        (3, 0, 0, 0, 0, 1): ([1, 7, 5, 9, 6, 8, 2, 3, 4, 10], [4, 10, 7, 6, 5, 8, 1, 2, 3, 9]),
    }
    for alpha, (v_lift, w_lift) in table.items():
        ul, vl, wl = block_lift(U, V, W, alpha)
        assert ul == U
        assert vl == Permutation(v_lift)
        assert wl == Permutation(w_lift)


def test_block_lift_lengths_and_bruhat_order():
    """(4.7): l(v~) = l(v) + l(sigma), l(w~) = l(w) + l(sigma) + |alpha|; v <= w implies v~ <= w~."""
    for alpha, (lv, lw) in {(2, 0, 0, 1, 0, 0): (16, 22), (1, 0, 0, 1, 0, 0): (11, 16), (3, 0, 0, 0, 0, 1): (19, 26)}.items():
        _, vl, wl = block_lift(U, V, W, alpha)
        assert vl.inv == lv and wl.inv == lw
        assert wl.inv - vl.inv == W.inv - V.inv + sum(alpha)
        assert vl.bruhat_leq(wl)
    s3 = Permutation.all_permutations(3)
    for v, w in itertools.product(s3, s3):
        for alpha in itertools.product(range(3), repeat=3):
            _, vl, wl = block_lift(Permutation([]), v, w, alpha)
            assert wl.inv - vl.inv == w.inv - v.inv + sum(alpha)
            if v.bruhat_leq(w):
                assert vl.bruhat_leq(wl)


def test_block_lift_zero_is_identity_and_pads():
    assert block_lift(U, V, W, (0,) * 6) == (U, V, W)
    # permutations shorter than alpha are padded with fixed points
    ul, vl, wl = block_lift([2, 1], [1, 3, 2], [], (1, 0, 2, 0))
    assert ul == Permutation([2, 1])
    assert vl == increasing_completion([1, 4, 3, 7], 7)
    assert wl == increasing_completion([2, 3, 6, 7], 7)


def test_block_lift_rejects_bad_input():
    with pytest.raises(ValueError):
        block_lift([2, 1], [1, 2], [1, 2], (1, -1))
    with pytest.raises(ValueError):
        block_lift([3, 1, 2], [1, 2], [1, 2], (1, 0))


def test_block_lift_schubert_coefficients_are_structure_constants():
    """Corollary 5.3 with v = w = id: [x^alpha] S_u(x) = c^{w~}_{u~,v~}, an ordinary structure constant."""
    sp = sympify_sympy
    n = 4
    gens = [sp(x[i]) for i in range(1, n + 1)]
    for u in Permutation.all_permutations(n):
        poly = sympy.Poly(sp(schubpoly(u, x, ZeroGeneratingSet())), *gens)
        for alpha in itertools.product(range(n), repeat=n):
            if sum(alpha) != u.inv:
                continue
            ul, vl, wl = block_lift(u, [], [], alpha)
            assert vl.bruhat_leq(wl)
            expected = poly.coeff_monomial(sympy.Mul(*[g**a for g, a in zip(gens, alpha)]))
            assert schubmult_py({ul: 1}, vl).get(wl, 0) == expected, (u, alpha)
