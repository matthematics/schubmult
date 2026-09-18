"""Tests for ``schubmult.mult.groth_quantum`` (single quantum Grothendieck kernel)."""

import itertools
import random

import pytest
import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y
from schubmult.mult.groth import grothmult_py
from schubmult.mult.groth_quantum import grothmult_q, grothmult_q_dict, grothmult_q_pieri
from schubmult.mult.groth_quantum_double import groth_elem_sym_poly_q, lm_quantize, quantum_pieri_chains
from schubmult.mult.quantum import schubmult_q
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars, grothendieck_poly

q = _vars.q_var
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


def _rat(rng):
    while True:
        r = sympy.Rational(rng.randint(-9, 9), rng.randint(1, 7))
        if r not in (0, 1, -1):
            return r


def test_quantum_pieri_chains_small():
    chains = quantum_pieri_chains(Permutation([2, 1]), 1)
    # empty chain, Bruhat cover (1,3), quantum edge (1,2) with weight q_1, and (1,3) then (1,2)
    assert chains[Permutation([2, 1])] == (0, (0, 0))
    assert chains[Permutation([3, 1, 2])] == (1, (0, 0))
    assert chains[Permutation([])] == (1, (1, 0))
    assert chains[Permutation([1, 3, 2])] == (2, (1, 0))
    assert len(chains) == 4


def test_quantum_pieri_chains_q0_is_classical_support():
    from schubmult.mult.groth_double import _top_block_support

    for u in S3:
        for k in (1, 2):
            classical = {w for w, (_, d) in quantum_pieri_chains(u, k).items() if not any(d)}
            assert classical == set(_top_block_support(u, k)) | {u}


def test_grothmult_q_known_values():
    res = grothmult_q({Permutation([2, 1]): S.One}, [2, 1])
    assert _same(res, {Permutation([3, 1, 2]): 1, Permutation([1, 3, 2]): beta * q[1], Permutation([]): q[1]})
    res = grothmult_q({Permutation([2, 3, 1]): S.One}, [3, 1, 2])
    assert _same(res, {Permutation([4, 2, 1, 3]): 1, Permutation([1, 2, 4, 3]): beta * q[1] * q[2], Permutation([]): q[1] * q[2]})


def test_grothmult_q_identity_and_dict():
    res = grothmult_q({Permutation([2, 1]): S.One}, [])
    assert _same(res, {Permutation([2, 1]): 1})
    d = grothmult_q_dict({Permutation([2, 1]): S.One}, {Permutation([2, 1]): 2, Permutation([]): 1})
    assert _same(d, {Permutation([3, 1, 2]): 2, Permutation([1, 3, 2]): 2 * beta * q[1], Permutation([]): 2 * q[1], Permutation([2, 1]): 1})


def test_grothmult_q_reduces_to_grothmult_py_at_q0():
    q0 = {sp(q[i]): 0 for i in range(1, 8)}
    for u in S3:
        for v in S3:
            assert _same(grothmult_q({u: S.One}, v), grothmult_py({u: S.One}, v), subs1=q0), (u, v)


def test_grothmult_q_reduces_to_schubmult_q_at_beta0():
    b0 = {sp(beta): 0}
    for u in S3:
        for v in S3:
            assert _same(grothmult_q({u: S.One}, v), schubmult_q({u: S.One}, v), subs1=b0), (u, v)


def test_grothmult_q_pieri_top_equals_product_with_s1():
    # Q(x_1) = G^q_{21}(x), so the top block for k = 1 is multiplication by G^q_{21}
    for u in S3:
        assert _same(grothmult_q_pieri({u: S.One}, 1, 1), grothmult_q({u: S.One}, [2, 1])), u


def _gq(w, bnum, qnum, cache):
    if w not in cache:
        y0 = {sp(y[i]): 0 for i in range(1, 12)} | {sp(beta): bnum}
        f = sp(grothendieck_poly(w, x, y, beta)).xreplace(y0)
        cache[w] = sp(lm_quantize(f, max(len(w), 2), x, bnum, qnum))
    return cache[w]


@pytest.mark.parametrize("seed", [1, 2])
def test_grothmult_q_pieri_polynomial_identity(seed):
    rng = random.Random(seed)
    bnum = _rat(rng)
    qnum = [None] + [_rat(rng) for _ in range(8)]
    spec = {sp(beta): bnum} | {sp(q[i]): qnum[i] for i in range(1, 9)}
    cache = {}
    for u in S3:
        for k in (1, 2):
            for p in range(1, k + 1):
                ep = sp(groth_elem_sym_poly_q(p, k, S.Zero, x, bnum, qnum))
                rhs = sum((sp(c).xreplace(spec) * _gq(w, bnum, qnum, cache) for w, c in grothmult_q_pieri({u: S.One}, p, k).items()), sympy.Integer(0))
                assert sympy.expand(ep * _gq(u, bnum, qnum, cache) - rhs) == 0, (u, p, k)


@pytest.mark.parametrize("seed", [3])
def test_grothmult_q_polynomial_identity(seed):
    rng = random.Random(seed)
    bnum = _rat(rng)
    qnum = [None] + [_rat(rng) for _ in range(8)]
    spec = {sp(beta): bnum} | {sp(q[i]): qnum[i] for i in range(1, 9)}
    cache = {}
    for u in S3:
        for v in S3:
            prod_dict = grothmult_q({u: S.One}, v)
            rhs = sum((sp(c).xreplace(spec) * _gq(w, bnum, qnum, cache) for w, c in prod_dict.items()), sympy.Integer(0))
            assert sympy.expand(_gq(u, bnum, qnum, cache) * _gq(v, bnum, qnum, cache) - rhs) == 0, (u, v)
