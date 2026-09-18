"""Tests for ``schubmult.mult.groth_quantum_double`` (quantum double Grothendieck kernel)."""

import itertools
import random
from math import comb

import pytest
import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y, z
from schubmult.mult.groth_double import grothmult_double, grothmult_double_pieri, grothmult_double_top
from schubmult.mult.groth_quantum import grothmult_q
from schubmult.mult.groth_quantum_double import (
    groth_elem_sym_poly_q,
    grothmult_q_double,
    grothmult_q_double_dict,
    grothmult_q_double_pieri,
    grothmult_q_double_top,
    lm_quantize,
    qgroth_poly,
    quantum_elem_sym,
    quantum_pieri_chains,
)
from schubmult.mult.quantum_double import mult_poly_q_double, schubmult_q_double
from schubmult.symbolic import S, Symbol, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars, elem_sym_poly_q, grothendieck_poly

q = _vars.q_var
zz = Symbol("zz")
S3 = [Permutation(list(p)) for p in itertools.permutations(range(1, 4))]
sp = sympify_sympy

Q0 = {sp(q[i]): 0 for i in range(1, 12)}
B0 = {sp(beta): 0}
FLIP_Y = {sp(y[i]): -sp(y[i]) for i in range(1, 12)}
FLIP_YZ = FLIP_Y | {sp(z[i]): -sp(z[i]) for i in range(1, 12)}


def _same(d1, d2, subs1=None, subs2=None):
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if subs2:
            b = b.xreplace(subs2)
        if sympy.cancel(a - b) != 0:
            return False
    return True


def _rat(rng):
    while True:
        r = sympy.Rational(rng.randint(-9, 9), rng.randint(1, 7))
        if r not in (0, 1, -1):
            return r


def _fgp_elem_sym(p, k, zvar):
    """beta = 0 limit of groth_elem_sym_poly_q(..., fgl=False): FGP quantization of e_p(x + z)."""
    xs = [x[i] for i in range(1, k + 1)]
    zeros = [S.Zero] * (k + 1)
    return sum((comb(k - j, p - j) * zvar ** (p - j) * elem_sym_poly_q(j, k, xs, zeros, q) for j in range(p + 1)), S.Zero)


# --- quantization -----------------------------------------------------------------------------


def test_quantum_elem_sym_degree_one():
    # F^2_1 = (1 - Q_1) X_1 + (1 - Q_2) X_2, Q_j = beta^2 q_j
    expected = (1 - beta**2 * q[1]) * (1 + beta * x[1]) + (1 - beta**2 * q[2]) * (1 + beta * x[2])
    assert sympy.expand(sp(quantum_elem_sym(1, 2, x, beta) - expected)) == 0


def test_lm_quantize_s1():
    g21 = grothendieck_poly(Permutation([2, 1]), x, y, beta)
    expected = x[1] + y[1] + beta * x[1] * y[1] - beta * q[1] * (1 + beta * x[1]) * (1 + beta * y[1])
    assert sympy.expand(sp(lm_quantize(g21, 2, x, beta) - expected)) == 0


def test_lm_quantize_is_linear_and_fixes_scalars():
    f = 3 * y[1] + beta * y[2] ** 2
    assert sympy.expand(sp(lm_quantize(f, 3, x, beta) - f)) == 0
    g = grothendieck_poly(Permutation([1, 3, 2]), x, y, beta)
    lhs = lm_quantize(y[1] * g + 2 * x[1], 3, x, beta)
    rhs = y[1] * lm_quantize(g, 3, x, beta) + 2 * lm_quantize(x[1], 3, x, beta)
    assert sympy.expand(sp(lhs - rhs)) == 0


def test_lm_quantize_stable_in_slots():
    g = grothendieck_poly(Permutation([2, 1]), x, y, beta)
    assert sympy.expand(sp(lm_quantize(g, 2, x, beta) - lm_quantize(g, 4, x, beta))) == 0


def test_lm_quantize_rejects_out_of_span():
    with pytest.raises(ValueError):
        lm_quantize(x[1] ** 2, 2, x, beta)


def test_qgroth_poly_matches_maeno_naito_sagaki_recursion():
    """At beta = -1 the LM quantization reproduces G^Q_w = pi^{(y)}_{w w0} G^Q_{w0} (MNS Part II)."""
    N = 3
    xs = [None] + [sympy.Symbol(f"x_{i}") for i in range(1, N + 1)]
    ys = [None] + [sympy.Symbol(f"y_{i}") for i in range(1, N + 1)]
    Qs = [None] + [sympy.Symbol(f"q_{i}") for i in range(1, N)] + [sympy.Integer(0)]

    def F(k, l):
        tot = 0
        for J in itertools.combinations(range(1, k + 1), l):
            term = 1
            for j in J:
                term *= 1 - xs[j]
                if j + 1 not in J:
                    term *= 1 - Qs[j]
            tot += term
        return sympy.expand(tot)

    def pi_y(i, f):
        s = f.subs({ys[i]: ys[i + 1], ys[i + 1]: ys[i]}, simultaneous=True)
        return sympy.cancel(f + (1 - ys[i]) * (f - s) / (ys[i] - ys[i + 1]))

    w0 = Permutation.w0(N)
    G0 = sympy.prod([sum((-1) ** l * (1 - ys[N - k]) ** l * F(k, l) for l in range(k + 1)) for k in range(1, N)])
    for w in S3:
        f = G0
        for i in reversed(list((w * w0).code_word)):
            f = pi_y(i, f)
        mine = sp(qgroth_poly(w, x, y, sympy.Integer(-1), q))
        assert sympy.expand(mine - f) == 0, w


def test_groth_elem_sym_poly_q_top_is_sum_of_elementaries():
    # Q(prod (x_i + z)) = sum_j z^{k-j} Q(e_j(x_1..x_k)); Q is linear over z
    k = 2
    top = groth_elem_sym_poly_q(k, k, zz, x, beta, fgl=False)
    rhs = sum((zz ** (k - j) * groth_elem_sym_poly_q(j, k, S.Zero, x, beta) for j in range(k + 1)), S.Zero)
    assert sympy.expand(sp(top - rhs)) == 0


# --- supports ---------------------------------------------------------------------------------


def test_quantum_pieri_chains_contains_empty_chain_and_is_deterministic():
    for u in S3:
        for k in (1, 2):
            chains = quantum_pieri_chains(u, k)
            assert chains[u][0] == 0
            assert not any(chains[u][1])
            assert chains is quantum_pieri_chains(u, k)


def test_quantum_pieri_chains_length_parity():
    # len - 2|D| = l(w) - l(u) along every quantum Bruhat chain
    for u in S3:
        for k in (1, 2):
            for w, (length, d) in quantum_pieri_chains(u, k).items():
                # each quantum edge (a, b) contributes b - a to |D| and drops l by 2(b - a) - 1
                assert length - 2 * sum(d) == w.inv - u.inv, (u, k, w)


# --- specializations of the rule ------------------------------------------------------------


def test_top_and_pieri_reduce_to_classical_at_q0():
    for u in S3:
        for k in (1, 2):
            assert _same(grothmult_q_double_top({u: S.One}, k, zz, y, beta), grothmult_double_top({u: S.One}, k, zz, y, beta), subs1=Q0), (u, k)
            for p in range(1, k + 1):
                assert _same(grothmult_q_double_pieri({u: S.One}, p, k, zz, y, beta), grothmult_double_pieri({u: S.One}, p, k, zz, None, y, beta), subs1=Q0), (u, p, k)


def test_top_and_pieri_reduce_to_quantum_schubert_at_beta0():
    # G^q_w(x; y)|_{beta=0} = S^q_w(x; -y)
    for u in S3:
        for k in (1, 2):
            top = _fgp_elem_sym(k, k, zz)
            assert _same(grothmult_q_double_top({u: S.One}, k, zz, y, S.Zero), mult_poly_q_double({u: S.One}, top, x, y, q), subs2=FLIP_Y), (u, k)
            for p in range(1, k + 1):
                ep = _fgp_elem_sym(p, k, zz)  # e_p(x (+) z) = e_p(x + z) at beta = 0
                assert _same(grothmult_q_double_pieri({u: S.One}, p, k, zz, y, S.Zero), mult_poly_q_double({u: S.One}, ep, x, y, q), subs2=FLIP_Y), (u, p, k)


def test_product_reduces_to_grothmult_double_at_q0():
    for u in S3:
        for v in S3:
            assert _same(grothmult_q_double({u: S.One}, v, y, z, beta), grothmult_double({u: S.One}, v, y, z, beta), subs1=Q0), (u, v)


def test_product_reduces_to_schubmult_q_double_at_beta0():
    for u in S3:
        for v in S3:
            assert _same(grothmult_q_double({u: S.One}, v, y, z, S.Zero), schubmult_q_double({u: S.One}, v, y, z, q), subs2=FLIP_YZ), (u, v)


def test_product_reduces_to_grothmult_q_at_yz0():
    yz0 = {sp(y[i]): 0 for i in range(1, 12)} | {sp(z[i]): 0 for i in range(1, 12)}
    for u in S3:
        for v in S3:
            assert _same(grothmult_q_double({u: S.One}, v, y, z, beta), grothmult_q({u: S.One}, v), subs1=yz0), (u, v)


def test_product_identity_and_dict():
    d = {Permutation([2, 1]): S.One}
    assert _same(grothmult_q_double(d, [], y, z, beta), d)
    single = grothmult_q_double(d, [2, 1], y, z, beta)
    combo = grothmult_q_double_dict(d, {Permutation([2, 1]): 3, Permutation([]): 1}, y, z, beta)
    expected = {w: 3 * c for w, c in single.items()}
    expected[Permutation([2, 1])] = expected.get(Permutation([2, 1]), 0) + 1
    assert _same(combo, expected)


def test_top_block_s1_known_value():
    # Q(x_1 + z) G^q_{21}(x; y).  The window value is u(1) = 2, so y_2 throughout: fixed on the
    # diagonal (z + (-)y_2), out for the Bruhat cover to 312, the quantum edge (1,2) to id (weight
    # q_1) and the length-2 chain (1,3),(1,2) to 132 (weight beta q_1).
    res = grothmult_q_double_top({Permutation([2, 1]): S.One}, 1, zz, y, beta)
    y2 = y[2]
    assert _same(
        res,
        {
            Permutation([2, 1]): (zz * (1 + beta * y2) - y2) / (1 + beta * y2),
            Permutation([3, 1, 2]): S.One / (1 + beta * y2),
            Permutation([]): q[1] / (1 + beta * y2),
            Permutation([1, 3, 2]): beta * q[1] / (1 + beta * y2),
        },
    )


# --- the polynomial identity ---------------------------------------------------------------


def _gq(w, var, spec, bnum, qnum, cache):
    key = (w, var.label)
    if key not in cache:
        f = sp(grothendieck_poly(w, x, var, beta)).xreplace(spec)
        cache[key] = sp(lm_quantize(f, max(len(w), 2), x, bnum, qnum))
    return cache[key]


@pytest.mark.parametrize("seed", [11, 12])
def test_polynomial_identity_random_specialization(seed):
    """G^q_u(x; y) G^q_v(x; z) = sum_w c^w_{uv} G^q_w(x; y) as polynomials in x, at random rational beta, y, z, q."""
    rng = random.Random(seed)
    bnum = _rat(rng)
    qnum = [None] + [_rat(rng) for _ in range(10)]
    spec = {sp(beta): bnum}
    spec |= {sp(y[i]): _rat(rng) for i in range(1, 10)}
    spec |= {sp(z[i]): _rat(rng) for i in range(1, 10)}
    spec |= {sp(q[i]): qnum[i] for i in range(1, 10)}
    cache = {}
    for u in S3:
        for v in S3:
            prod_dict = grothmult_q_double({u: S.One}, v, y, z, beta, q)
            rhs = sum((sp(c).xreplace(spec) * _gq(w, y, spec, bnum, qnum, cache) for w, c in prod_dict.items()), sympy.Integer(0))
            lhs = _gq(u, y, spec, bnum, qnum, cache) * _gq(v, z, spec, bnum, qnum, cache)
            assert sympy.expand(lhs - rhs) == 0, (u, v)


def test_top_block_polynomial_identity_random_specialization():
    rng = random.Random(5)
    bnum = _rat(rng)
    qnum = [None] + [_rat(rng) for _ in range(10)]
    znum = _rat(rng)
    spec = {sp(beta): bnum, sp(zz): znum}
    spec |= {sp(y[i]): _rat(rng) for i in range(1, 10)}
    spec |= {sp(q[i]): qnum[i] for i in range(1, 10)}
    cache = {}
    for u in S3:
        for k in (1, 2):
            top = sp(groth_elem_sym_poly_q(k, k, znum, x, bnum, qnum, fgl=False))
            rhs = sum((sp(c).xreplace(spec) * _gq(w, y, spec, bnum, qnum, cache) for w, c in grothmult_q_double_top({u: S.One}, k, zz, y, beta).items()), sympy.Integer(0))
            assert sympy.expand(top * _gq(u, y, spec, bnum, qnum, cache) - rhs) == 0, (u, k)
