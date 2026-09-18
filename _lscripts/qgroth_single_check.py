"""Checks for the single quantum Grothendieck kernel ``schubmult.mult.groth_quantum``.

    python _lscripts/qgroth_single_check.py [n]

1. q = 0 reduces to ``grothmult_py``; beta = 0 reduces to ``schubmult_q``.
2. ``grothmult_q_double`` at y = z = 0 agrees with ``grothmult_q``.
3. Polynomial identity ``G^q_u G^q_v = sum c^w G^q_w`` with ``G^q = lm_quantize(G)`` (symbolic).
4. beta = -1: ``grothmult_q_pieri`` reproduces the Naito--Sagaki quantum K Pieri theorem.
"""

import itertools
import sys
import time

import sympy

from schubmult import Permutation
from schubmult.abc import beta, x, y, z
from schubmult.combinatorics.permutation import uncode
from schubmult.mult.groth import grothmult_py
from schubmult.mult.groth_quantum import grothmult_q, grothmult_q_pieri
from schubmult.mult.groth_quantum_double import grothmult_q_double, lm_quantize
from schubmult.mult.quantum import schubmult_q
from schubmult.symbolic import S, sympify_sympy
from schubmult.symbolic.poly.schub_poly import _vars, grothendieck_poly

n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
q = _vars.q_var
perms = [Permutation(list(p)) for p in itertools.permutations(range(1, n + 1))]
t0 = time.time()
sp = sympify_sympy
q0 = {sp(q[i]): 0 for i in range(1, 3 * n)}
yz0 = {sp(y[i]): 0 for i in range(1, 3 * n)} | {sp(z[i]): 0 for i in range(1, 3 * n)}


def dict_diff(d1, d2, subs1=None, subs2=None):
    bad = []
    for w in set(d1) | set(d2):
        a, b = sp(d1.get(w, 0)), sp(d2.get(w, 0))
        if subs1:
            a = a.xreplace(subs1)
        if subs2:
            b = b.xreplace(subs2)
        if sympy.expand(a - b) != 0:
            bad.append(w)
    return bad


fails_q0 = fails_b0 = 0
for u in perms:
    for v in perms:
        res = grothmult_q({u: S.One}, v, beta, q)
        if dict_diff(res, grothmult_py({u: S.One}, v, beta), subs1=q0):
            fails_q0 += 1
        if dict_diff(res, schubmult_q({u: S.One}, v), subs1={sp(beta): 0}):
            fails_b0 += 1
print(f"S_{n}: q=0 vs grothmult_py {fails_q0} failures; beta=0 vs schubmult_q {fails_b0} failures   ({time.time() - t0:.0f}s)", flush=True)

fails = 0
for u in perms:
    for v in perms:
        if dict_diff(grothmult_q_double({u: S.One}, v, y, z, beta, q), grothmult_q({u: S.One}, v, beta, q), subs1=yz0):
            fails += 1
print(f"S_{n}: grothmult_q_double at y=z=0 vs grothmult_q: {fails} failures   ({time.time() - t0:.0f}s)", flush=True)

cache = {}


def gq(w):
    if w not in cache:
        f = sp(grothendieck_poly(w, x, y, beta)).xreplace(yz0)
        cache[w] = sp(lm_quantize(f, max(len(w), 2), x, beta, q))
    return cache[w]


fails = 0
for u in perms:
    for v in perms:
        res = grothmult_q({u: S.One}, v, beta, q)
        rhs = sum((sp(c) * gq(w) for w, c in res.items()), sympy.Integer(0))
        if sympy.expand(gq(u) * gq(v) - rhs) != 0:
            fails += 1
            if fails <= 3:
                print("  polynomial identity fails", list(u), list(v))
print(f"S_{n}: polynomial identity G^q_u G^q_v: {len(perms) ** 2} products, {fails} failures   ({time.time() - t0:.0f}s)", flush=True)


# Naito--Sagaki (arXiv:2211.01578, Thm 2.13) at beta = -1, Q_j = q_j
def ns_prec(l1, l2):
    (a, b), (c, d) = l1, l2
    return b > d or (b == d and a < c)


def ns_pieri_chains(u, k, N):
    labels = [(a, b) for b in range(N, k, -1) for a in range(1, k + 1)]
    out = []

    def walk(w, path, qmon):
        out.append((tuple(path), w, qmon))
        for lab in labels:
            a, b = lab
            if lab in path or (path and b > path[-1][1]):
                continue
            if len(path) >= 2 and any(path[t][0] == path[-1][0] for t in range(len(path) - 1)) and not ns_prec(path[-1], lab):
                continue
            w2 = w.swap(a - 1, b - 1)
            if w2.inv == w.inv + 1:
                walk(w2, [*path, lab], qmon)
            elif w2.inv == w.inv - (2 * (b - a) - 1):
                walk(w2, [*path, lab], qmon * sympy.prod([sp(q[j]) for j in range(a, b)]))

    walk(u, [], sympy.Integer(1))
    return out


def ns_markings(path, p):
    r = len(path)
    forced = set()
    t = 0
    while t < r and all(path[i][1] == path[0][1] for i in range(t + 1)) and all(path[i][0] > path[i + 1][0] for i in range(t)):
        forced.add(t)
        t += 1
    count = 0
    for Mk in itertools.combinations(range(r), p):
        Ms = set(Mk)
        if not forced <= Ms:
            continue
        ok = True
        for s in range(r):
            if s in Ms:
                if any(path[t_][0] == path[s][0] for t_ in range(s)):
                    ok = False
                    break
            elif s < r - 1 and not ns_prec(path[s], path[s + 1]):
                ok = False
                break
        if ok:
            count += 1
    return count


# G^q_{c[k,p]} (beta = -1) in terms of Q(e_j(x_1..x_k)): G_{c[k,p]}(x) = sum_j (-1)^{j-p} binom(j-1, p-1) e_j  (since
# G_{c[k,p]} = e_p(x) + (-1) e_{p+1}... check numerically instead: expand G_{c[k,p]} in e_j via lm_quantize-free algebra)
def groth_cyclic_in_elem(k, p):
    """Coefficients c_j with G_{c[k,p]}(x)|_{beta=-1} = sum_j c_j e_j(x_1..x_k), by solving on symmetric polynomials."""
    xs = [sp(x[i]) for i in range(1, k + 1)]
    cpk = uncode([0] * (k - p) + [1] * p)
    f = sympy.expand(sp(grothendieck_poly(cpk, x, y, beta)).xreplace(yz0).xreplace({sp(beta): -1}))
    es = [sum(sympy.prod(J) for J in itertools.combinations(xs, j)) for j in range(k + 1)]
    cs = sympy.symbols(f"c0:{k + 1}")
    sol = sympy.solve(sympy.Poly(f - sum(c * e for c, e in zip(cs, es)), *xs).coeffs(), cs, dict=True)
    assert len(sol) == 1, (k, p, sol)
    return [sol[0].get(c, 0) for c in cs]


fails = tot = 0
N = n
for u in perms:
    for k in range(1, N):
        chains = ns_pieri_chains(u, k, N + k)
        cyc = {p: groth_cyclic_in_elem(k, p) for p in range(1, k + 1)}
        pieri = {j: grothmult_q_pieri({u: S.One}, j, k, sympy.Integer(-1), q) for j in range(0, k + 1)}
        for p in range(1, k + 1):
            pred = {}
            for path, w, qmon in chains:
                mk = ns_markings(path, p)
                if mk:
                    pred[w] = pred.get(w, 0) + (-1) ** (len(path) - p) * mk * qmon
            got = {}
            for j, cj in enumerate(cyc[p]):
                if cj != 0:
                    for w, c in pieri[j].items():
                        got[w] = got.get(w, 0) + cj * sp(c)
            tot += 1
            if any(sympy.expand(got.get(w, 0) - pred.get(w, 0)) != 0 for w in set(got) | set(pred)):
                fails += 1
                if fails <= 3:
                    print("  NS mismatch", list(u), k, p, {tuple(w): c for w, c in got.items()}, {tuple(w): c for w, c in pred.items()})
print(f"S_{n}: Naito--Sagaki theorem via grothmult_q_pieri at beta=-1: {tot} products, {fails} mismatches   ({time.time() - t0:.0f}s)")
