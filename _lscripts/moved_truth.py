"""HONEST ground-truth verifier for the moved case (u NOT->_k w, u(k)!=w(k)).

Builds the ACTUAL recursion RHS F as a polynomial (factorial elementary syms, like
moved_inspect2.py), groups by q-monomial, and checks expand(sum)==0 for EVERY q-monomial,
over ALL moved cases and all admissible p. No combinatorial up/down counting proxy.

This is the only trustworthy F=0 test. Reports total groups checked and any nonzero.
"""
import sys
from collections import defaultdict
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

y = GeneratingSet("y"); z = GeneratingSet("z"); q_gs = GeneratingSet("q")
def q_ab(a, b): return prod([q_gs[s] for s in range(a, b)])
def val(perm, pos): return perm[pos-1]

def enumerate_pieri(u, k, N):
    u = Permutation(u); results = {u: S.One}
    stack = [(u, frozenset(), N+1, S.One, u.inv)]; seen = set()
    while stack:
        perm, used_a, last_b, qw, clen = stack.pop(); key = (perm, used_a, last_b)
        if key in seen: continue
        seen.add(key)
        for a in range(1, k+1):
            if a in used_a: continue
            for b in range(k+1, N+1):
                if b > last_b: continue
                nperm = perm.swap(a-1, b-1); d = nperm.inv - clen
                if d == 1: nqw = qw
                elif d == -2*(b-a)+1: nqw = qw*q_ab(a, b)
                else: continue
                results.setdefault(nperm, nqw)
                stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def Pset(u, w, m): return sorted(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))

def cinfo(u, w, m, a, Rm):
    """c_m(a;u,w) = q_m(u,w) * E_{a-n_m, m-n_m}(y_P; z), nonzero only if u->_m w. Returns
    (qmon, poly) or None."""
    if m < 0: return None
    w = Permutation(w)
    if w not in Rm: return None
    qm = Rm[w]; Pm = Pset(u, w, m); nm = m-len(Pm); deg = a-nm; nv = m-nm
    if deg < 0 or deg > nv: return None
    if nv == 0:
        return (qm, S.One) if deg == 0 else None
    yv = [y[v] for v in Pm]; nc = nv+1-deg; zv = [z[i] for i in range(1, nc+1)]
    return (qm, expand_func(FactorialElemSym(deg, nv, yv, zv)))

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    moved = 0; groups_checked = 0; nonzero = 0; examples = []
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N)
            Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue           # moved case: u NOT->_k w
                moved += 1
                for p in range(1, k+1):
                    groups = defaultdict(lambda: S.Zero)
                    def add(qmon, poly):
                        groups[str(expand(qmon))] += poly
                    ci = cinfo(u, w, k-1, p-1, Rkm1)
                    if ci: add(ci[0], (y[val(w, k)]-z[k-p+1])*ci[1])   # DIAG
                    ci = cinfo(u, w, k-1, p, Rkm1)
                    if ci: add(ci[0], ci[1])                            # STRAIGHT (REC)
                    ci = cinfo(u, w, k-2, p-2, Rkm2)
                    if ci: add(q_gs[k-1]*ci[0], ci[1])                  # QUANT
                    for pp in range(1, N+1):
                        if pp == k: continue
                        vv = w.swap(k-1, pp-1); dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1: extra = S.One
                        elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
                        else: continue
                        sgn = 1 if pp > k else -1
                        ci = cinfo(u, vv, k-1, p-1, Rkm1)
                        if ci: add(extra*ci[0], sgn*ci[1])              # UP / DOWN Monk
                    for kk, poly in groups.items():
                        groups_checked += 1
                        if expand(poly) != 0:
                            nonzero += 1
                            if len(examples) < 20:
                                examples.append((tuple(u), tuple(w), k, p, kk, expand(poly)))
    print(f"n={n}  moved cases={moved}  q-monomial groups checked={groups_checked}")
    print(f"  NONZERO groups (real F_M != 0) = {nonzero}")
    for ex in examples:
        print(f"   u={ex[0]} w={ex[1]} k={ex[2]} p={ex[3]} qmon={ex[4]} sum={ex[5]}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
