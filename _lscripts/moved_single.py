"""Exact F for a single (u,w,k) across all p, fully expand-based, listing every term."""
import sys
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, expand_func, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import FactorialElemSym
y = GeneratingSet("y"); z = GeneratingSet("z"); q_gs = GeneratingSet("q")
def q_ab(a, b): return prod([q_gs[s] for s in range(a, b)])
def val(p, i): return p[i-1]
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
                results.setdefault(nperm, nqw); stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results
def Pset(u, w, m): return sorted(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def cpoly(u, w, m, a, Rm):
    if m < 0 or w not in Rm: return None
    qm = Rm[w]; Pm = Pset(u, w, m); nm = m-len(Pm); deg = a-nm; nv = m-nm
    if deg < 0 or deg > nv: return None
    if nv == 0: return (qm, S.One) if deg == 0 else None
    yv = [y[v] for v in Pm]; nc = nv+1-deg; zv = [z[i] for i in range(1, nc+1)]
    return (qm, expand_func(FactorialElemSym(deg, nv, yv, zv)))
def main():
    u = Permutation(list(map(int, sys.argv[1].split(","))))
    w = Permutation(list(map(int, sys.argv[2].split(","))))
    k = int(sys.argv[3]); N = int(sys.argv[4]) if len(sys.argv) > 4 else 2*len(u)+4
    Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N); Rkm2 = enumerate_pieri(u, k-2, N)
    print(f"u={tuple(u)} w={tuple(w)} k={k} N={N}")
    print(f"  u->_k w ? {w in Rk}   u->_{{k-1}} w ? {w in Rkm1}   u->_{{k-2}} w ? {w in Rkm2}")
    print(f"  u(k)={val(u,k)} w(k)={val(w,k)}")
    for p in range(1, k+1):
        terms = []
        ci = cpoly(u, w, k-1, p-1, Rkm1)
        if ci: terms.append(("DIAG", (y[val(w,k)]-z[k-p+1])*ci[0]*ci[1]))
        ci = cpoly(u, w, k-1, p, Rkm1)
        if ci: terms.append(("STRAIGHT", ci[0]*ci[1]))
        ci = cpoly(u, w, k-2, p-2, Rkm2)
        if ci: terms.append(("QUANT", q_gs[k-1]*ci[0]*ci[1]))
        for pp in range(1, N+1):
            if pp == k: continue
            vv = w.swap(k-1, pp-1); dl = w.inv-vv.inv; lo, hi = min(pp, k), max(pp, k)
            if dl == 1: extra = S.One
            elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
            else: continue
            ci = cpoly(u, vv, k-1, p-1, Rkm1)
            if not ci: continue
            sgn = 1 if pp > k else -1
            terms.append((("UP" if pp > k else "DOWN")+f"(b={pp})", sgn*extra*ci[0]*ci[1]))
        F = expand(sum(t for _, t in terms))
        print(f"  p={p}: F={F}")
        for nm, t in terms:
            print(f"      {nm}: {expand(t)}")
if __name__ == "__main__":
    main()
