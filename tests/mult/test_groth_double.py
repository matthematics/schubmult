"""Tests for ``schubmult.mult.groth_double`` (double beta-Grothendieck Monk/Pieri formula)."""

import itertools

import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y, z
from schubmult.mult.double import schubmult_double
from schubmult.mult.groth_double import (
    _tilde_elem_sym_frac,
    _top_block_support,
    dgroth_phantom_expansion,
    double_groth_times_double_schub,
    epsilon_chain,
    grothmult_double,
    grothmult_double_pieri,
    grothmult_double_top,
    groth_elem_sym_func,
    monk_chain,
    mult_poly_groth_double,
    single_variable_groth,
)
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import grothendieck_poly, schubpoly
from schubmult.symbolic.poly.variables import ZeroGeneratingSet
from schubmult.utils.test_utils import vanishes

zero = ZeroGeneratingSet()
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy
FLIP_YZ = {sp(y[i]): -sp(y[i]) for i in range(1, 12)} | {sp(z[i]): -sp(z[i]) for i in range(1, 12)}


def _groth(v, var2, var3):
    return sp(grothendieck_poly(v, var2, var3, beta))


def _same(d1, d2, subs1=None, subs2=None):
    # coefficients are rational functions in beta/y/z: test the differences exactly at random points
    diffs = []
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if subs2:
            b = b.xreplace(subs2)
        diffs.append(a - b)
    return vanishes(diffs)


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
            assert vanishes([_groth(u, x, y) * _groth(v, x, z) - rhs]), (u, v)


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


def test_double_groth_times_double_schub_polynomial_identity():
    # G_u(x; y) S_v(x; z) = sum_w c_w G_w(x; y) as rational functions; S_v(x; z) = schubpoly(v, x, z).
    from schubmult.symbolic.poly.schub_poly import schubpoly

    S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
    for u in S4:
        for v in S3:
            d = double_groth_times_double_schub(u, v, y, z, beta)
            lhs = _groth(u, x, y) * sp(schubpoly(v, x, z))
            rhs = sum((sp(c) * _groth(w, x, y) for w, c in d.items()), sympy.Integer(0))
            assert vanishes([lhs - rhs]), (u, v)


def test_double_groth_times_double_schub_defaults_and_identity_perm():
    u = Permutation([2, 1, 3])
    assert _same(double_groth_times_double_schub(u, Permutation([])), {u: S.One})
    assert _same(double_groth_times_double_schub(u, Permutation([1, 3, 2])), double_groth_times_double_schub(u, Permutation([1, 3, 2]), y, z, beta))


def _conjugate(th):
    return [sum(1 for part in th if part >= i) for i in range(1, th[0] + 1)]


def _formal_inverse_x():
    return [None] + [-x[i] / (1 + beta * x[i]) for i in range(1, 8)]


def test_dgroth_phantom_expansion_identity():
    # G_v(x; z) = prod (1 + beta x_i)^{lambda'_i} sum_sigma c_sigma S_sigma((-)x; z), lambda = theta(v^{-1}).
    S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
    xm = _formal_inverse_x()
    for v in S4:
        th = [t for t in (~v).theta() if t]
        exp = dgroth_phantom_expansion(v, z, beta)
        if not th:
            assert exp == {v: S.One}
            continue
        prefactor = sympy.prod([(1 + beta * x[i]) ** e for i, e in enumerate(_conjugate(th), start=1)])
        rhs = prefactor * sum((sp(c) * sp(schubpoly(sg, xm, z)) for sg, c in exp.items()), sympy.Integer(0))
        assert vanishes([_groth(v, x, z) - rhs]), v


def test_dgroth_phantom_expansion_strict_theta_identity():
    # The same identity over the strict-theta shape used by the quantum kernel.
    S4 = [Permutation(list(p)) for p in itertools.permutations(range(1, 5))]
    xm = _formal_inverse_x()
    for v in S4:
        th = [t for t in (~v).strict_theta() if t]
        if not th:
            continue
        exp = dgroth_phantom_expansion(v, z, beta, th)
        prefactor = sympy.prod([(1 + beta * x[i]) ** e for i, e in enumerate(_conjugate(th), start=1)])
        rhs = prefactor * sum((sp(c) * sp(schubpoly(sg, xm, z)) for sg, c in exp.items()), sympy.Integer(0))
        assert vanishes([_groth(v, x, z) - rhs]), v


def test_dgroth_phantom_expansion_rejects_bad_shape():
    import pytest

    with pytest.raises(ValueError):
        dgroth_phantom_expansion(Permutation([2, 3, 1]), z, beta, [1])


def test_tilde_elem_sym_frac_matches_exact_multiplication():
    # Coefficient of G_{u2} in prod_{j<=k}(1 + beta x_j) E_{p,k}((-)x; z_sel) G_{u1}, against the exact fold.
    from schubmult.mult.groth_double import _frac_to_expr

    def layer_factor(p, k, zs):
        # (-1)^p sum_{|I|=p} prod_{j not in I}(1 + beta x_j) prod_m (x_{i_m} (+) z_{i_m - m + 1})
        total = S.Zero
        for chosen in itertools.combinations(range(1, k + 1), p):
            term = S.One
            for j in range(1, k + 1):
                if j not in chosen:
                    term *= S.One + beta * x[j]
            for m, i in enumerate(chosen, start=1):
                term *= x[i] + zs[i - m] + beta * x[i] * zs[i - m]
            total += term
        return (-1) ** p * total

    from schubmult import uncode
    from schubmult.symbolic.poly.schub_poly import call_zvars
    from schubmult.utils.schub_lib import compute_vpathdicts

    # genuine v-path steps (v1 -> v2, vdiff) of a few v, with their layer data (k, i)
    steps = set()
    for v in (Permutation([2, 3, 1]), Permutation([3, 1, 2]), Permutation([1, 4, 2, 3]), Permutation([2, 4, 1, 3])):
        th = [t for t in (~v).theta() if t]
        for index, layer in enumerate(compute_vpathdicts(tuple(th), v * uncode(th))):
            for v1, moves in layer.items():
                for v2, vdiff, _ in moves:
                    steps.add((th[index], index + 1, v1, v2, vdiff))
    for u1 in S3:
        for k, i, v1, v2, vdiff in sorted(steps, key=str):
            zidx = call_zvars(v1, v2, k, i)[: vdiff + 1]
            exact = mult_poly_groth_double({u1: S.One}, layer_factor(k - vdiff, k, [z[c] for c in zidx]), x, y, beta)
            closed = {u2: _frac_to_expr(_tilde_elem_sym_frac(k, i, u1, u2, v1, v2, vdiff, y, z, beta), y, beta) for u2 in _top_block_support(u1, k) | {u1}}
            assert _same(exact, closed), (u1, k, i, v1, v2, vdiff)
