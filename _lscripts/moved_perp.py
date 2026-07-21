"""Affirmative per-p verification of the proof-supporting facts (expand-based).

For each moved case and each admissible p, partition the recursion terms by q-monomial.
For every PURE reflection group (no DIAG/STRAIGHT/QUANT at this p) verify:
  (P1) all p-admissible members share the SAME (P_{k-1}, n_{k-1})  => identical E-coefficient;
  (P2) #up == #down among p-admissible members;
  (P3) expand(group sum) == 0.
Also confirm groups containing DIAG/STRAIGHT contain NO up-Monk term (Case i is a pure
telescope), and QUANT groups have exactly one Monk partner. Report all fail counts.
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
    if m < 0: return None
    w = Permutation(w)
    if w not in Rm: return None
    qm = Rm[w]; Pm = Pset(u, w, m); nm = m-len(Pm); deg = a-nm; nv = m-nm
    if deg < 0 or deg > nv: return None
    if nv == 0:
        return (qm, S.One, (tuple(Pm), nm)) if deg == 0 else None
    yv = [y[v] for v in Pm]; nc = nv+1-deg; zv = [z[i] for i in range(1, nc+1)]
    return (qm, expand_func(FactorialElemSym(deg, nv, yv, zv)), (tuple(Pm), nm))

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    pure_groups = 0; failP1 = failP2 = failP3 = 0
    special_with_up = 0; quant_partner_bad = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N); Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rk) | set(Rkm1) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue
                for p in range(1, k+1):
                    groups = defaultdict(list)
                    kind = defaultdict(set)   # qkey -> set of {'diag','straight','quant'}
                    ci = cinfo(u, w, k-1, p-1, Rkm1)
                    if ci:
                        qk = str(expand(ci[0])); groups[qk].append(("DIAG", 1, (y[val(w, k)]-z[k-p+1])*ci[1], ci[2])); kind[qk].add("ds")
                    ci = cinfo(u, w, k-1, p, Rkm1)
                    if ci:
                        qk = str(expand(ci[0])); groups[qk].append(("STRAIGHT", 1, ci[1], ci[2])); kind[qk].add("ds")
                    ci = cinfo(u, w, k-2, p-2, Rkm2)
                    if ci:
                        qk = str(expand(q_gs[k-1]*ci[0])); groups[qk].append(("QUANT", 1, ci[1], ci[2])); kind[qk].add("quant")
                    for pp in range(1, N+1):
                        if pp == k: continue
                        vv = w.swap(k-1, pp-1); dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1: extra = S.One
                        elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
                        else: continue
                        sgn = 1 if pp > k else -1
                        ci = cinfo(u, vv, k-1, p-1, Rkm1)
                        if ci:
                            qk = str(expand(extra*ci[0]))
                            groups[qk].append((("UP" if pp > k else "DOWN"), sgn, ci[1], ci[2]))
                    for qk, mem in groups.items():
                        ups = [m for m in mem if m[0] == "UP"]
                        downs = [m for m in mem if m[0] == "DOWN"]
                        if "ds" in kind[qk]:
                            if ups: special_with_up += 1     # Case (i) telescope should have no up
                            continue
                        if "quant" in kind[qk]:
                            monk = ups + downs
                            if not (len(ups) == 0 and len(downs) == 1):
                                quant_partner_bad += 1
                            continue
                        # pure reflection group
                        pure_groups += 1
                        pn = set(m[3] for m in mem)
                        if len(pn) != 1: failP1 += 1
                        if len(ups) != len(downs): failP2 += 1
                        if expand(sum(s*poly for _, s, poly, _ in mem)) != 0: failP3 += 1
    print(f"n={n} pure reflection groups (per-p) = {pure_groups}")
    print(f"  (P1) members not sharing (P,n) : fails = {failP1}")
    print(f"  (P2) #up != #down              : fails = {failP2}")
    print(f"  (P3) expand(sum) != 0          : fails = {failP3}")
    print(f"  diag/straight group containing an UP term : {special_with_up}")
    print(f"  quant group without unique DOWN partner   : {quant_partner_bad}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
