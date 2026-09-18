"""Equivariant quantum K oracle for Fl_N from the Maeno--Naito--Sagaki presentation, and the
equivariant quantum top-block experiment.

QK_T(Fl_N) = K[x_1..x_N] / I^Q,  K = coefficient field in y_1..y_N, Q_1..Q_{N-1},
    I^Q = ( F^N_l(x) - e_l(1/(1-y_1), ..., 1/(1-y_N)) : l = 1..N ),
    F^k_l(x) = sum_{J in [k], |J|=l} prod_{j in J, j+1 not in J} (1 - Q_j) prod_{j in J} (1 - x_j),  1 - Q_N := 1
(MNS Part I, arXiv:2302.09485, Thm 6.1; e^{-eps_j} = 1 - y_j).  Schubert classes are the
Lenart--Maeno quantum double Grothendieck polynomials (MNS Part II, arXiv:2305.17685):
    G^Q_{w0}(x,y) = prod_{k=1}^{N-1} sum_{l=0}^{k} (-1)^l (1 - y_{N-k})^l F^k_l(x),
    G^Q_w = pi^{(y)}_{w w0} G^Q_{w0},   pi^{(y)}_i f = f + (1 - y_i)(f - s_i f)/(y_i - y_{i+1}).
Everything is at beta = -1; the beta-generic statement follows by homogeneity (deg x = deg y = 1,
deg beta = -1, deg Q = 2).

Usage:  python qk_equivariant_oracle.py N [sym|num] [trials]
  sym: all parameters symbolic (feasible for N = 3);
  num: y, Q, z specialized to random rationals, several trials (default; feasible for N = 4).

Checks: Q = 0 reproduces DGx (beta = -1); y = 0 reproduces the Naito--Sagaki Pieri theorem.
Experiment: coefficients of G^Q_w(x;y) in  Q(prod_{i<=k}(x_i + z)) * G^Q_u(x;y), the quantized top
block sum_l (-1)^l (1+z)^{k-l} F^k_l(x), against
    (-1)^(len-m) Q^D  prod_Fix (z - y_a/(1-y_a))  prod_Left (1+z)  prod_Out (1-y_a)^(-1)
over the Naito--Sagaki marked k-Pieri chains in the quantum Bruhat graph.
"""

import itertools
import random
import sys
import time

import sympy
from sympy.polys.matrices import DomainMatrix

from schubmult import DGx, Permutation
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import sympify_sympy

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3
MODE = sys.argv[2] if len(sys.argv) > 2 else ("sym" if N <= 3 else "num")
TRIALS = int(sys.argv[3]) if len(sys.argv) > 3 else 3

x = [None] + list(sympy.symbols(f"x1:{N + 1}"))
y = [None] + list(sympy.symbols(f"y1:{N + 1}"))
Q = [None] + list(sympy.symbols(f"Q1:{N}")) + [sympy.Integer(0)]  # Q_N = 0 so that 1 - Q_N = 1
z = sympy.Symbol("z")

perms = sorted((Permutation(list(p)) for p in itertools.permutations(range(1, N + 1))), key=lambda w: (w.inv, tuple(w)))
idx = {w: i for i, w in enumerate(perms)}
M = len(perms)
w0 = Permutation.w0(N)
t0 = time.time()


def F(k, l):
    tot = 0
    for J in itertools.combinations(range(1, k + 1), l):
        Js = set(J)
        term = 1
        for j in J:
            term *= 1 - x[j]
            if j + 1 not in Js:
                term *= 1 - Q[j]
        tot += term
    return sympy.expand(tot)


# ---- quantum double Grothendieck polynomials (symbolic) -------------------------------------
def pi_y(i, f):
    s = f.subs({y[i]: y[i + 1], y[i + 1]: y[i]}, simultaneous=True)
    return sympy.cancel(f + (1 - y[i]) * (f - s) / (y[i] - y[i + 1]))


G0 = sympy.prod([sum((-1) ** l * (1 - y[N - k]) ** l * F(k, l) for l in range(k + 1)) for k in range(1, N)])
GQ_sym = {}
for w in perms:
    f = G0
    for i in reversed(list((w * w0).code_word)):  # w w0 = s_{i_1}...s_{i_l}; innermost first
        f = pi_y(i, f)
    GQ_sym[w] = f
print(f"quantum double Grothendieck polynomials in {time.time() - t0:.1f}s", flush=True)

bad = 0
sub_beta = {sympify_sympy(DGx([2, 1]).ring._beta): -1}
subQ0 = {Q[j]: 0 for j in range(1, N)}
for w in perms:
    ref = sympify_sympy(DGx(w).expand()).subs(sub_beta)
    ref = ref.subs({sympy.Symbol(f"x_{j}"): x[j] for j in range(1, N + 1)} | {sympy.Symbol(f"y_{j}"): y[j] for j in range(1, N + 1)})
    if sympy.cancel(GQ_sym[w].subs(subQ0) - ref) != 0:
        bad += 1
print(f"G^Q_w at Q=0 vs DGx(beta=-1): {bad} mismatches", flush=True)


# ---- the quotient ring for a given specialization of the parameters ---------------------------
class Oracle:
    def __init__(self, spec):
        """spec: dict of parameter substitutions (may be empty for the symbolic run)."""
        self.spec = spec
        free = [p for p in y[1:] + Q[1:N] + [z] if p not in spec]
        self.K = sympy.QQ.frac_field(*free) if free else sympy.QQ
        yv = [None] + [y[j].subs(spec) for j in range(1, N + 1)]
        gens = [F(N, l).subs(spec) - sum(sympy.prod([1 / (1 - yv[j]) for j in J]) for J in itertools.combinations(range(1, N + 1), l)) for l in range(1, N + 1)]
        self.GB = sympy.groebner(gens, *x[1:], order="grevlex", domain=self.K)
        self.GQ = {w: sympy.expand(GQ_sym[w].subs(spec)) for w in perms}
        nfs = {w: sympy.Poly(self.nf(self.GQ[w]), *x[1:], domain=self.K) for w in perms}
        monos = sorted({m for p in nfs.values() for m in p.monoms()})
        assert len(monos) == M, (len(monos), M)
        self.mono_idx = {m: i for i, m in enumerate(monos)}
        rows = [[self.K.zero] * M for _ in range(M)]
        for w in perms:
            for m, c in nfs[w].terms():
                rows[self.mono_idx[m]][idx[w]] = self.K.convert(c)
        self.B_inv = DomainMatrix(rows, (M, M), self.K).inv()

    def nf(self, expr):
        return self.GB.reduce(sympy.expand(expr))[1]

    def coords(self, expr):
        p = sympy.Poly(self.nf(expr), *x[1:], domain=self.K)
        col = [[self.K.zero] for _ in range(M)]
        for m, c in p.terms():
            col[self.mono_idx[m]][0] = self.K.convert(c)
        c = (self.B_inv * DomainMatrix(col, (M, 1), self.K)).rep.to_list()
        return {w: self.K.to_sympy(c[i][0]) for i, w in enumerate(perms) if c[i][0] != self.K.zero}

    def product(self, u, mult):
        return self.coords(self.GQ[u] * sympy.sympify(mult).subs(self.spec))


# ---- Naito--Sagaki k-Pieri chains and markings ----------------------------------------------
def ns_prec(l1, l2):
    (a, b), (c, d) = l1, l2
    return b > d or (b == d and a < c)


def ns_pieri_chains(u, k):
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
                walk(w2, path + [lab], qmon)
            elif w2.inv == w.inv - (2 * (b - a) - 1):
                walk(w2, path + [lab], qmon * sympy.prod([Q[j] for j in range(a, b)]))

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
                if any(path[u_][0] == path[s][0] for u_ in range(s)):
                    ok = False
                    break
            elif s < r - 1 and not ns_prec(path[s], path[s + 1]):
                ok = False
                break
        if ok:
            count += 1
    return count


def fate(u, w, k):
    window = [w[i] for i in range(k)]
    fix, left, out = [], [], []
    for i in range(k):
        v = u[i]
        if window[i] == v:
            fix.append(v)
        elif v in window and window.index(v) < i:
            left.append(v)
        else:
            out.append(v)
    return fix, left, out


def quantum_chain_data(u, k):
    data = {}
    for path, w, qmon in ns_pieri_chains(u, k):
        if any(ns_markings(path, p) for p in range(0, k + 1)):
            data.setdefault(w, set()).add((len(path), qmon))
    return data


def check_ns(oracle):
    bad_ns = tot = 0
    for u in perms:
        for k in range(1, N):
            chains = ns_pieri_chains(u, k)
            for p in range(0, k + 1):
                cpk = uncode([0] * (k - p) + [1] * p)
                pred = {}
                for path, w, qmon in chains:
                    mk = ns_markings(path, p)
                    if mk:
                        pred[w] = pred.get(w, 0) + (-1) ** (len(path) - p) * mk * qmon.subs(oracle.spec)
                got = oracle.product(u, GQ_sym[cpk])
                tot += 1
                if any(sympy.simplify(got.get(w, 0) - pred.get(w, 0)) != 0 for w in set(got) | set(pred)):
                    bad_ns += 1
                    if bad_ns <= 3:
                        print("  NS mismatch", u, k, p, {tuple(w): c for w, c in got.items()}, {tuple(w): c for w, c in pred.items()})
    return bad_ns, tot


def check_top_block(oracle, verbose):
    agree = disagree = 0
    for u in perms:
        for k in range(1, N):
            mult = sum((-1) ** l * (1 + z) ** (k - l) * F(k, l) for l in range(k + 1))
            got = oracle.product(u, mult)
            data = quantum_chain_data(u, k)
            for w in sorted(set(got) | set(data), key=lambda p: (p.inv, tuple(p))):
                fix, left, out = fate(u, w, k)
                m = len(left) + len(out)
                preds = set()
                for ln, qmon in data.get(w, set()):
                    pr = (-1) ** (ln - m) * qmon * sympy.prod([z - y[a] / (1 - y[a]) for a in fix]) * (1 + z) ** len(left) * sympy.prod([1 / (1 - y[a]) for a in out])
                    preds.add(sympy.cancel(pr.subs(oracle.spec)))
                gotc = sympy.cancel(got.get(w, 0))
                if len(preds) == 1 and sympy.cancel(gotc - next(iter(preds))) == 0:
                    agree += 1
                else:
                    disagree += 1
                    if verbose or disagree <= 10:
                        print(f"  u={list(u)} k={k} w={list(w)}: got {sympy.factor(gotc)}  predicted {[sympy.factor(p) for p in preds]}  fate fix={fix} left={left} out={out}", flush=True)
    return agree, disagree


if MODE == "sym":
    oracle = Oracle({})
    print(f"symbolic oracle built in {time.time() - t0:.0f}s", flush=True)
    b, t = check_ns(Oracle({y[j]: 0 for j in range(1, N + 1)}))
    print(f"Naito--Sagaki theorem at y=0: {t} products, {b} mismatches   ({time.time() - t0:.0f}s)", flush=True)
    a, d = check_top_block(oracle, verbose=True)
    print(f"\nequivariant quantum top block (symbolic): {a} coefficients agree, {d} do not   ({time.time() - t0:.0f}s)")
else:
    random.seed(11)
    tot_a = tot_d = 0
    for trial in range(TRIALS):
        def rat():
            while True:
                r = sympy.Rational(random.randint(-19, 19), random.randint(2, 11))
                if r not in (0, 1):
                    return r
        spec_ns = {y[j]: 0 for j in range(1, N + 1)} | {Q[j]: rat() for j in range(1, N)} | {z: rat()}
        if trial == 0:
            b, t = check_ns(Oracle(spec_ns))
            print(f"Naito--Sagaki theorem at y=0, random Q: {t} products, {b} mismatches   ({time.time() - t0:.0f}s)", flush=True)
        spec = {y[j]: rat() for j in range(1, N + 1)} | {Q[j]: rat() for j in range(1, N)} | {z: rat()}
        oracle = Oracle(spec)
        a, d = check_top_block(oracle, verbose=False)
        tot_a += a
        tot_d += d
        print(f"trial {trial + 1}: spec {dict((str(k), v) for k, v in spec.items())}: {a} agree, {d} disagree   ({time.time() - t0:.0f}s)", flush=True)
    print(f"\nequivariant quantum top block, {TRIALS} random specializations: {tot_a} agree, {tot_d} disagree")
