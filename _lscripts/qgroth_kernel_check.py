"""Checks for the quantum double Grothendieck kernel ``schubmult.mult.groth_quantum_double``.

    python _lscripts/qgroth_kernel_check.py [n]

1. q = 0: ``grothmult_q_double_top`` / ``grothmult_q_double_pieri`` reduce to the classical kernels.
2. beta = 0: they reduce to the quantum double Schubert kernels (``mult_poly_q_double``, ``schubmult_q_double``).
3. Polynomial identity: ``G^q_u(x; y) G^q_v(x; z) = sum_w c^w_{uv} G^q_w(x; y)`` in ``Z[beta, q, y, z][x]``.
4. beta = -1: ``qgroth_poly`` agrees with the Maeno--Naito--Sagaki quantum double Grothendieck recursion.
"""

import itertools
import sys
import time

import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y, z
from schubmult.mult.groth_double import grothmult_double, grothmult_double_pieri, grothmult_double_top
from schubmult.mult.groth_quantum_double import grothmult_q_double, grothmult_q_double_pieri, grothmult_q_double_top, qgroth_poly
from schubmult.mult.quantum_double import mult_poly_q_double, schubmult_q_double
from schubmult.symbolic import S, Symbol, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars, elem_sym_poly_q
from math import comb

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
q = _vars.q_var
zz = Symbol("zz")
perms = [Permutation(list(p)) for p in itertools.permutations(range(1, n + 1))]
t0 = time.time()


def sp(e):
    return sympify_sympy(e)


def dict_diff(d1, d2, subs=None, subs2=None):
    bad = []
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs:
            a, b = a.subs(subs), b.subs(subs)
        if subs2:
            b = b.subs(subs2, simultaneous=True)
        if sympy.cancel(a - b) != 0:
            bad.append(w)
    return bad


def fgp_elem_sym(p, k, zvar, fgl):
    """beta = 0 limit of groth_elem_sym_poly_q: FGP quantization of e_p(x + z)."""
    xs = [x[i] for i in range(1, k + 1)]
    zeros = [S.Zero] * (k + 1)
    return sum((comb(k - j, p - j) * zvar ** (p - j) * elem_sym_poly_q(j, k, xs, zeros, q) for j in range(p + 1)), S.Zero)


q0 = {sp(q[i]): 0 for i in range(1, 2 * n + 2)}
flip_y = {sp(y[i]): -sp(y[i]) for i in range(1, 2 * n + 2)}
flip_yz = flip_y | {sp(z[i]): -sp(z[i]) for i in range(1, 2 * n + 2)}

# 1. q = 0
fails = 0
for u in perms:
    for k in range(1, n):
        if dict_diff(grothmult_q_double_top({u: S.One}, k, zz, y, beta), grothmult_double_top({u: S.One}, k, zz, y, beta), q0):
            fails += 1
        for p in range(1, k + 1):
            if dict_diff(grothmult_q_double_pieri({u: S.One}, p, k, zz, y, beta), grothmult_double_pieri({u: S.One}, p, k, zz, None, y, beta), q0):
                fails += 1
print(f"q=0 top/pieri vs classical kernels on S_{n}: {fails} failures   ({time.time() - t0:.0f}s)", flush=True)

# 2. beta = 0  (G^q_w(x; y) at beta = 0 is S^q_w(x; -y), hence the sign flips)
fails = 0
for u in perms:
    for k in range(1, n):
        top = fgp_elem_sym(k, k, zz, False)
        if dict_diff(grothmult_q_double_top({u: S.One}, k, zz, y, S.Zero), mult_poly_q_double({u: S.One}, top, x, y, q), subs2=flip_y):
            fails += 1
        for p in range(1, k + 1):
            ep = fgp_elem_sym(p, k, zz, True)
            if dict_diff(grothmult_q_double_pieri({u: S.One}, p, k, zz, y, S.Zero), mult_poly_q_double({u: S.One}, ep, x, y, q), subs2=flip_y):
                fails += 1
print(f"beta=0 top/pieri vs mult_poly_q_double on S_{n}: {fails} failures   ({time.time() - t0:.0f}s)", flush=True)

fails = 0
for u in perms:
    for v in perms:
        if dict_diff(grothmult_q_double({u: S.One}, v, y, z, S.Zero), schubmult_q_double({u: S.One}, v, y, z, q), subs2=flip_yz):
            fails += 1
            if fails <= 3:
                print("  beta=0 mismatch", list(u), list(v))
print(f"beta=0 full product vs schubmult_q_double on S_{n}: {fails} failures   ({time.time() - t0:.0f}s)", flush=True)

# 1b. q = 0 full product
fails = 0
for u in perms:
    for v in perms:
        if dict_diff(grothmult_q_double({u: S.One}, v, y, z, beta), grothmult_double({u: S.One}, v, y, z, beta), q0):
            fails += 1
print(f"q=0 full product vs grothmult_double on S_{n}: {fails} failures   ({time.time() - t0:.0f}s)", flush=True)

# 4. beta = -1 vs the MNS recursion (symbolic; the recursion itself is cheap for n <= 4)
N = n
xs = [None] + [sympy.Symbol(f"x_{i}") for i in range(1, N + 1)]
ys = [None] + [sympy.Symbol(f"y_{i}") for i in range(1, N + 1)]
Qs = [None] + [sympy.Symbol(f"q_{i}") for i in range(1, N)] + [sympy.Integer(0)]


def F(k, l):
    tot = 0
    for J in itertools.combinations(range(1, k + 1), l):
        Js = set(J)
        term = 1
        for j in J:
            term *= 1 - xs[j]
            if j + 1 not in Js:
                term *= 1 - Qs[j]
        tot += term
    return sympy.expand(tot)


def pi_y(i, f):
    s = f.subs({ys[i]: ys[i + 1], ys[i + 1]: ys[i]}, simultaneous=True)
    return sympy.cancel(f + (1 - ys[i]) * (f - s) / (ys[i] - ys[i + 1]))


w0 = Permutation.w0(N)
G0 = sympy.prod([sum((-1) ** l * (1 - ys[N - k]) ** l * F(k, l) for l in range(k + 1)) for k in range(1, N)])
fails = 0
for w in perms:
    f = G0
    for i in reversed(list((w * w0).code_word)):
        f = pi_y(i, f)
    mine = sp(qgroth_poly(w, x, y, sympy.Integer(-1), q))
    if sympy.expand(mine - f) != 0:
        fails += 1
        if fails <= 3:
            print("  MNS mismatch", list(w), sympy.factor(mine), sympy.factor(f))
print(f"beta=-1 qgroth_poly vs Maeno--Naito--Sagaki G^Q_w on S_{n}: {fails} mismatches   ({time.time() - t0:.0f}s)", flush=True)

# 3. polynomial identity, under random rational specializations of beta, y, z, q (x symbolic)
import random  # noqa: E402

from schubmult.symbolic.poly.schub_poly import grothendieck_poly  # noqa: E402
from schubmult.mult.groth_quantum_double import lm_quantize  # noqa: E402

random.seed(7)
TRIALS = int(sys.argv[2]) if len(sys.argv) > 2 else 2


def rat():
    while True:
        r = sympy.Rational(random.randint(-9, 9), random.randint(1, 7))
        if r not in (0, 1, -1):
            return r


for trial in range(TRIALS):
    spec = {sp(beta): rat()}
    spec |= {sp(y[i]): rat() for i in range(1, 3 * n)}
    spec |= {sp(z[i]): rat() for i in range(1, 3 * n)}
    spec |= {sp(q[i]): rat() for i in range(1, 3 * n)}
    bnum = spec[sp(beta)]
    qnum = [None] + [spec[sp(q[i])] for i in range(1, 3 * n)]
    cache = {}

    def gq(w, var):
        key = (w, var.label)
        if key not in cache:
            f = sp(grothendieck_poly(w, x, var, beta)).xreplace(spec)
            cache[key] = sp(lm_quantize(f, max(len(w), 2), x, bnum, qnum))
        return cache[key]

    fails = 0
    count = 0
    for u in perms:
        for v in perms:
            prod_dict = grothmult_q_double({u: S.One}, v, y, z, beta, q)
            rhs = sum((sp(c).xreplace(spec) * gq(w, y) for w, c in prod_dict.items()), sympy.Integer(0))
            count += 1
            if sympy.expand(gq(u, y) * gq(v, z) - rhs) != 0:
                fails += 1
                if fails <= 3:
                    print("  polynomial identity fails", list(u), list(v))
    print(f"trial {trial + 1}: polynomial identity G^q_u(x;y) G^q_v(x;z) on S_{n}: {count} products, {fails} failures   ({time.time() - t0:.0f}s)", flush=True)
