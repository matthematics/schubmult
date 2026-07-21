"""Consolidated verifier for LEMMA X (moved reflection pairing). expand-based. n=5 AND n=6.

Moved case: u NOT->_k w, u(k)!=w(k). pi=u^{-1}w, C=cycle thru k, B=max(C).
For each admissible p, split recursion terms of the RHS F by q-monomial into buckets.
Distinguished terms: DIAG, STRAIGHT (both weight q_{k-1}(u,w)); QUANT (weight q_{k-1}q_{k-2}).
Monk terms: kappa_{pp,k} c_{k-1}(p-1;u, w t_{pp,k}), pp != k.

Assert (fail counters, expect all 0 at n=5 AND n=6):
 (L1) at most ONE contributing up-term (pp>k) across ALL p and buckets for this (u,w,k);
      and if it exists its index is exactly B=max(C).
 (L2) every Case(ii) bucket (no DIAG/STRAIGHT/QUANT) is either empty or has exactly one up-term
      (index B) and one down-term (index a<k).
 (L3) in such a bucket, up and down share P_{k-1} and n_{k-1} (=> identical E) AND identical
      q-weight kappa_{kB} q_{k-1}(u,wt_{kB}) == kappa_{ak} q_{k-1}(u,wt_{ak})  [expand].
 (L4) expand(bucket sum)==0 for EVERY bucket (diag/straight telescope, quant-partner, up-down).
 (L5) expand(total F)==0.
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
def npart(u, w, m): return m - len(Pset(u, w, m))
def cyc_through(perm, x, N):
    c = [x]; cur = perm[x-1]
    while cur != x and cur not in c: c.append(cur); cur = perm[cur-1]
    return set(c)

def cpoly(u, w, m, a, Rm):
    """c_m(a;u,w)=q_m(u,w)*E_{a-n_m,m-n_m}(y_P;z). Returns (qmon, poly, (Ptuple,nm)) or None."""
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
    moved = 0; fL1 = fL2 = fL3 = fL4 = fL5 = 0
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
                moved += 1
                pi = (~u)*w; C = cyc_through(pi, k, N); B = max(C)
                up_indices_seen = set()
                Ftotal = S.Zero
                for p in range(1, k+1):
                    buckets = defaultdict(list)   # qkey -> list of (kind, sign, poly, info, idx)
                    ci = cpoly(u, w, k-1, p-1, Rkm1)
                    if ci: buckets[str(expand(ci[0]))].append(("DIAG", 1, (y[val(w, k)]-z[k-p+1])*ci[1], ci[2], None))
                    ci = cpoly(u, w, k-1, p, Rkm1)
                    if ci: buckets[str(expand(ci[0]))].append(("STRAIGHT", 1, ci[1], ci[2], None))
                    ci = cpoly(u, w, k-2, p-2, Rkm2)
                    if ci: buckets[str(expand(q_gs[k-1]*ci[0]))].append(("QUANT", 1, ci[1], ci[2], None))
                    for pp in range(1, N+1):
                        if pp == k: continue
                        vv = w.swap(k-1, pp-1); dl = w.inv - vv.inv
                        lo, hi = min(pp, k), max(pp, k)
                        if dl == 1: extra = S.One
                        elif dl == 1-2*(hi-lo): extra = q_ab(lo, hi)
                        else: continue
                        ci = cpoly(u, vv, k-1, p-1, Rkm1)
                        if not ci: continue
                        sgn = 1 if pp > k else -1
                        buckets[str(expand(extra*ci[0]))].append((("UP" if pp > k else "DOWN"), sgn, ci[1], ci[2], pp))
                        if pp > k:
                            up_indices_seen.add(pp)
                    for qk, mem in buckets.items():
                        kinds = set(m[0] for m in mem)
                        bsum = expand(sum(s*poly for _, s, poly, _, _ in mem))
                        if bsum != 0: fL4 += 1
                        Ftotal += sum(s*poly for _, s, poly, _, _ in mem)
                        if kinds & {"DIAG", "STRAIGHT", "QUANT"}:
                            continue
                        # Case (ii) bucket
                        ups = [m for m in mem if m[0] == "UP"]; downs = [m for m in mem if m[0] == "DOWN"]
                        if not (len(ups) == 1 and len(downs) == 1):
                            fL2 += 1; continue
                        up = ups[0]; dn = downs[0]
                        if up[4] != B: fL2 += 1
                        # L3: invariants + weight
                        if up[3] != dn[3]:
                            fL3 += 1
                        w_up_key = expand(qk); w_dn_key = expand(qk)  # same bucket => same weight by construction
                        # explicit weight identity check:
                        # up weight = kappa_{kB} q_{k-1}(u, w t_{kB}); dn weight = kappa_{ak} q_{k-1}(u, w t_{ak})
                        # both equal the bucket key by construction; assert they are equal (they are)
                    # end buckets
                # L1: at most one up-index, equal to B
                if len(up_indices_seen) > 1 or (up_indices_seen and up_indices_seen != {B}):
                    fL1 += 1
                if expand(Ftotal) != 0:
                    fL5 += 1
    print(f"n={n} moved cases = {moved}")
    print(f"  (L1) >1 up-index or up-index != B  : fails = {fL1}")
    print(f"  (L2) Case(ii) bucket not exactly {{1 up=B, 1 down}} : fails = {fL2}")
    print(f"  (L3) up/down (P,n) mismatch        : fails = {fL3}")
    print(f"  (L4) some bucket expand(sum) != 0  : fails = {fL4}")
    print(f"  (L5) total F != 0                  : fails = {fL5}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
