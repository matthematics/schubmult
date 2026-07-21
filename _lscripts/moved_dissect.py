"""Per-p dissection of pure-reflection groups (no diag/straight/quant at that p) to find the
REAL cancellation mechanism. Uses expand on actual factorial-elementary polynomials.

For each moved case, fixed p, and each q-monomial group that contains ONLY Monk (up/down)
terms (no DIAG/STRAIGHT/QUANT), print the p-admissible members (label, sign, (deg,nv,P),
poly) and the expanded sum. Focus on groups where the count of up-admissible != down-admissible
to expose whether cancellation is (a) pairing, (b) recurrence among several, or (c) E=0.
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
        return (qm, S.One, (deg, nv, tuple(Pm))) if deg == 0 else None
    yv = [y[v] for v in Pm]; nc = nv+1-deg; zv = [z[i] for i in range(1, nc+1)]
    return (qm, expand_func(FactorialElemSym(deg, nv, yv, zv)), (deg, nv, tuple(Pm)))

def run(n, want=8):
    N = 2*n; perms = list(Permutation.all_permutations(n)); shown = 0
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
                    groups = defaultdict(list)   # qkey -> list of (label, sign, poly, info)
                    has_special = defaultdict(bool)
                    ci = cinfo(u, w, k-1, p-1, Rkm1)
                    if ci: groups[str(expand(ci[0]))].append(("DIAG", 1, (y[val(w, k)]-z[k-p+1])*ci[1], ci[2])); has_special[str(expand(ci[0]))] = True
                    ci = cinfo(u, w, k-1, p, Rkm1)
                    if ci: groups[str(expand(ci[0]))].append(("STRAIGHT", 1, ci[1], ci[2])); has_special[str(expand(ci[0]))] = True
                    ci = cinfo(u, w, k-2, p-2, Rkm2)
                    if ci:
                        qk = str(expand(q_gs[k-1]*ci[0])); groups[qk].append(("QUANT", 1, ci[1], ci[2])); has_special[qk] = True
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
                            lab = f"{'UP' if pp>k else 'DOWN'}(r={pp})"
                            groups[qk].append((lab, sgn, ci[1], ci[2]))
                    for qk, mem in groups.items():
                        if has_special[qk]: continue   # pure reflection groups only
                        ups = [m for m in mem if m[0].startswith("UP")]
                        downs = [m for m in mem if m[0].startswith("DOWN")]
                        if len(ups) == len(downs): continue   # naive-balanced; skip
                        tot = expand(sum(s*poly for _, s, poly, _ in mem))
                        shown += 1
                        print(f"\n=== u={tuple(u)} w={tuple(w)} k={k} p={p} qmon={qk} "
                              f"#up={len(ups)} #down={len(downs)} sum={'0' if tot==0 else tot} ===")
                        for lab, s, poly, info in mem:
                            print(f"    {lab:12s} sign={s:+d} (deg,nv,P)={info}  poly={expand(poly)}")
                        if shown >= want: return
    print("done shown", shown)

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
