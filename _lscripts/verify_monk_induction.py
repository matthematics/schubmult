"""Test the inductive Monk argument for  (*)  N(u,mu_v) <= Lambda(mu_v) N(u,v).

Setup: fix v; mu = mu_v = dom(v^{-1})^{-1} (dominant hull, LEFT weak order).
  Lam_mu = |[e,mu]_L|,  Lam_v = |[e,v]_L|.
  N(u,w) = |U_{w'<=_L w} supp(S_u S_{w'})|.
  R(u,v) = N(u,mu)/N(u,v).   Want R <= Lam_mu.   Base R(e,v)=Lam_mu/Lam_v.

User's argument (induction on ell(u) via Monk / Bruhat covers u0 <. u):
  (B1) numerator Lipschitz:   N(u,mu) <= (n-1) N(u0,mu)  ?
  (B2) denominator growth:    N(u,v) >= N(u0,v)          (monotone) ?
  (claim1) N(s_k,mu) <= (n-1) N(1,mu) = (n-1) Lam_mu     ?
  (claim2) N(s_k,v)  >= (n-1)  "usually"                 ?
  (STEP-strong) exists right descent s (u0=us <. u) with R(u,v) <= R(u0,v) ?
                (if always, proves R <= Lam_mu/Lam_v, the STRONGER P1p -- known to FAIL)
"""
import sys
from itertools import combinations
from schubmult import Permutation
from schubmult.rings.schubert.schubert_ring import Sx


def left_downset(perm, n):
    start = tuple(perm[i] for i in range(n))
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        pos = [0] * (n + 2)
        for idx, val in enumerate(a):
            pos[val] = idx
        for i in range(1, n):
            if pos[i] > pos[i + 1]:
                b = list(a)
                b[pos[i]] = i + 1
                b[pos[i + 1]] = i
                b = tuple(b)
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
    return seen


_supp = {}


def supp(u, w):
    key = (u, w)
    r = _supp.get(key)
    if r is None:
        r = frozenset((Sx(Permutation(list(u))) * Sx(Permutation(list(w)))).keys())
        _supp[key] = r
    return r


def length(a, n):
    return sum(1 for i in range(n) for j in range(i + 1, n) if a[i] > a[j])


def run(n):
    perms = [tuple(p[i] for i in range(n)) for p in Permutation.all_permutations(n)]
    # precompute downsets & N for every (u,target-perm)
    ds = {p: left_downset(Permutation(list(p)), n) for p in perms}
    Ncache = {}

    def N(u, w):
        key = (u, w)
        r = Ncache.get(key)
        if r is None:
            acc = set()
            for wp in ds[w]:
                acc |= supp(u, tuple(wp))
            r = len(acc)
            Ncache[key] = r
        return r

    # mu_v for each v
    mu_of = {}
    Lam = {p: len(ds[p]) for p in perms}
    for p in perms:
        vi = ~Permutation(list(p))
        mu = ~(vi.minimal_dominant_above())
        mu_of[p] = tuple(mu[i] for i in range(n))

    simple = []
    for k in range(n - 1):
        a = list(range(1, n + 1))
        a[k], a[k + 1] = a[k + 1], a[k]
        simple.append(tuple(a))

    e = tuple(range(1, n + 1))

    # --- claim1 / claim2 over simple reflections and all v ---
    max_c1 = 0.0     # N(s_k,mu)/Lam_mu   (want <= n-1)
    min_c2 = 10**9   # N(s_k,v)           (want >= n-1 "usually")
    c2_exceptions = 0
    c2_total = 0
    exc_examples = []
    # NEW: restrict to ell(v) >= n  (exclude "short" v)
    c2_exc_long = 0
    min_c2_long = 10**9
    exc_long_examples = []
    for v in perms:
        mu = mu_of[v]
        Lm = Lam[mu]
        lv = length(v, n)
        for s in simple:
            max_c1 = max(max_c1, N(s, mu) / Lm)
            nv = N(s, v)
            c2_total += 1
            if nv < n - 1:
                c2_exceptions += 1
                if len(exc_examples) < 6:
                    exc_examples.append((s, v, nv, N(s, mu), Lm))
                if lv >= n:
                    c2_exc_long += 1
                    if len(exc_long_examples) < 8:
                        exc_long_examples.append((s, v, lv, nv, N(s, mu), Lm))
            min_c2 = min(min_c2, nv)
            if lv >= n:
                min_c2_long = min(min_c2_long, nv)

    # --- (B1) numerator Lipschitz & (B2) denom monotone over Bruhat covers ---
    # covers: u0 <. u  iff  u = u0 * t_{ab}, ell(u)=ell(u0)+1
    covers = []  # (u0, u)
    permset = set(perms)
    for u in perms:
        lu = length(u, n)
        for a, b in combinations(range(n), 2):
            w = list(u)
            w[a], w[b] = w[b], w[a]
            w = tuple(w)
            if length(w, n) == lu - 1:
                covers.append((w, u))  # w <. u

    maxB1 = 0.0
    B2_ok = True
    step_strong_all = True   # for every u!=e exists cover with R(u)<=R(u0)
    # organize covers by top u
    from collections import defaultdict
    cov_by_top = defaultdict(list)
    for (u0, u) in covers:
        cov_by_top[u].append(u0)

    worstR = (0.0, None)
    for v in perms:
        mu = mu_of[v]
        Lm = Lam[mu]
        Rof = {}
        for u in perms:
            nu_mu = N(u, mu)
            nu_v = N(u, v)
            Rof[u] = nu_mu / nu_v
            if nu_mu / Lm > worstR[0]:
                worstR = (nu_mu / (Lm * nu_v), (u, v, mu, nu_mu, nu_v, Lm, Lam[v]))
        for u in perms:
            if u == e:
                continue
            # B1, B2 across its covers
            good_exists = False
            for u0 in cov_by_top[u]:
                if N(u0, mu) > 0:
                    maxB1 = max(maxB1, N(u, mu) / N(u0, mu))
                if N(u, v) < N(u0, v):
                    B2_ok = False
                if Rof[u] <= Rof[u0] + 1e-9:
                    good_exists = True
            if not good_exists:
                step_strong_all = False

    print(f"===== n={n} =====")
    print(f"  (claim1) max N(s_k,mu)/Lam_mu = {max_c1:.4f}   (<= n-1={n-1}? {'YES' if max_c1<=n-1+1e-9 else 'NO'})")
    print(f"  (claim2) min N(s_k,v) = {min_c2};  exceptions N(s_k,v)<n-1: {c2_exceptions}/{c2_total}")
    print(f"  (claim2 | ell(v)>=n) min N(s_k,v) = {min_c2_long};  exceptions: {c2_exc_long}")
    for (s, v, lv, nv, nmu, Lm) in exc_long_examples:
        print(f"       [ell(v)={lv}] s_k={s} v={v}: N(s,v)={nv}  N(s,mu)={nmu}  Lam_mu={Lm}  "
              f"(*)ratio={nmu/(Lm*nv):.3f}")
    print(f"  (B1) max N(u,mu)/N(u0,mu) over covers = {maxB1:.4f}   (<= n-1={n-1}? {'YES' if maxB1<=n-1+1e-9 else 'NO'})")
    print(f"  (B2) N(.,v) monotone up Bruhat covers: {B2_ok}")
    print(f"  (STEP-strong) every u has cover with R(u)<=R(u0): {step_strong_all}  (False => needs slack)")
    r, info = worstR
    print(f"  (*) max N(u,mu)/(Lam_mu N(u,v)) = {r:.4f}   {info}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
