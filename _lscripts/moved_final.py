"""FINAL structural verification for the moved case u NOT->_k w with u(k)!=w(k).

Claims (expect 0 fails, n=5 AND n=6):
 (A) u NOT->_{k-1} w  and  u NOT->_{k-2} w.
     [=> diagonal, straight, and quantum terms of lemma:vanishing all vanish here.]
 (B) The cycle C of pi=u^{-1}w through k has exactly two indices >= k: k and B=max C.
 (C) contributing up-neighbours = { t_{kB} } at most (any contributing up has b == B);
     #contrib-up in {0,1}, #contrib-down in {0,1}, and #contrib-up == #contrib-down.
 (D) when both == 1 (up v_up=w t_{kB}, down v_dn=w t_{rk}, r<k):
       P_{k-1}(u,v_up) == P_{k-1}(u,v_dn) == P_{k-1}(u,w),
       n_{k-1} equal,
       kappa_{kB} q_{k-1}(u,v_up) == kappa_{rk} q_{k-1}(u,v_dn)   (expand-equal),
     hence the up Monk term (sign +) and down Monk term (sign -) cancel:
       expand( kappa_up * c_{k-1}(p-1;u,v_up) )  ==  expand( kappa_dn * c_{k-1}(p-1;u,v_dn) )
     for every p (checked via full factorial-elementary E, ground truth).
"""
import sys
from schubmult import *  # noqa
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, prod
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
                results.setdefault(nperm, nqw); stack.append((nperm, used_a | {a}, b, nqw, nperm.inv))
    return results

def Pset(u, w, m): return frozenset(val(u, i) for i in range(1, m+1) if val(u, i) == val(w, i))
def npart(u, w, m): return m - len(Pset(u, w, m))
def cyc_through(perm, x):
    c = [x]; cur = perm[x-1]
    while cur != x and cur not in c: c.append(cur); cur = perm[cur-1]
    return c

def Efac(a, m, P):
    # E_{a,m}(y_P; z): factorial elementary, |P| variables. 0 if a<0 or a>m.
    if a < 0 or a > m: return S.Zero
    if a == 0: return S.One
    Pl = sorted(P)
    return FactorialElemSym(a, m, [y[p] for p in Pl], z)

def cco(u, v, m, a, Rm):
    # c_m(a;u,v) = q_m E_{a-n, m-n}(y_P; z), 0 if u NOT->_m v
    if v not in Rm: return S.Zero
    n = npart(u, v, m); P = Pset(u, v, m)
    return Rm[v]*Efac(a-n, m-n, P)

def run(n):
    N = 2*n; perms = list(Permutation.all_permutations(n))
    moved = 0; fA = fB = fC = fD = fCancel = 0
    for u in perms:
        u = Permutation(u)
        for k in range(2, n):
            Rk = enumerate_pieri(u, k, N); Rkm1 = enumerate_pieri(u, k-1, N); Rkm2 = enumerate_pieri(u, k-2, N)
            cands = set(Rkm1) | set(Rk) | set(Rkm2)
            for v0 in list(Rkm1):
                for pp in range(1, N+1):
                    if pp != k: cands.add(v0.swap(k-1, pp-1))
            for w in cands:
                w = Permutation(w)
                if val(u, k) == val(w, k): continue
                if w in Rk: continue  # only u NOT->_k w
                moved += 1
                pi = (~u)*w; C = cyc_through(pi, k); B = max(C)
                # (A)
                if (w in Rkm1) or (w in Rkm2): fA += 1
                # (B)
                geqk = sorted(x for x in C if x >= k)
                if geqk != sorted({k, B}): fB += 1
                # contributing ups/downs
                ups = []
                for b in range(k+1, N+1):
                    vv = w.swap(k-1, b-1); dl = w.inv - vv.inv
                    if dl == 1: kap = S.One
                    elif dl == 1-2*(b-k): kap = q_ab(k, b)
                    else: continue
                    if vv in Rkm1: ups.append((b, kap, vv))
                downs = []
                for a in range(1, k):
                    vv = w.swap(a-1, k-1); dl = w.inv - vv.inv
                    if dl == 1: kap = S.One
                    elif dl == 1-2*(k-a): kap = q_ab(a, k)
                    else: continue
                    if vv in Rkm1: downs.append((a, kap, vv))
                # (C)
                if not (len(ups) <= 1 and len(downs) <= 1 and len(ups) == len(downs)
                        and all(b == B for b, _, _ in ups)):
                    fC += 1
                # (D)
                if len(ups) == 1 and len(downs) == 1:
                    b, kap_up, v_up = ups[0]; a, kap_dn, v_dn = downs[0]
                    Pu = Pset(u, v_up, k-1); Pd = Pset(u, v_dn, k-1); Pw = Pset(u, w, k-1)
                    if not (Pu == Pd == Pw and npart(u, v_up, k-1) == npart(u, v_dn, k-1)):
                        fD += 1
                    if expand(kap_up*Rkm1[v_up]) != expand(kap_dn*Rkm1[v_dn]):
                        fD += 1
                    # full cancellation per p
                    for p in range(1, k+1):
                        up_term = kap_up*cco(u, v_up, k-1, p-1, Rkm1)
                        dn_term = kap_dn*cco(u, v_dn, k-1, p-1, Rkm1)
                        if expand(up_term - dn_term) != S.Zero:
                            fCancel += 1; break
    print(f"n={n}: moved cases (u NOT->_k w, u(k)!=w(k)) = {moved}")
    print(f"  (A) u NOT->_{{k-1}}w and NOT->_{{k-2}}w : fails = {fA}")
    print(f"  (B) C has exactly two indices>=k {{k,B}} : fails = {fB}")
    print(f"  (C) #up,#down in{{0,1}}, equal, up at B  : fails = {fC}")
    print(f"  (D) P/n/weight match up<->down          : fails = {fD}")
    print(f"  (cancel) up Monk term == down Monk term  : fails = {fCancel}")

if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
