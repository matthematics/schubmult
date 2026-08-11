"""Verify the sharp clean lemma behind the length threshold:

  (B+)  N(u,v) >= Lambda_v := |[e,v]_L|   for ALL u,v   (strengthens fact B; here for u=s_k and all u)
  and   Lambda_v >= ell(v)+1              (chain in weak order)

Consequences for the denominator claim:
  N(s_k,v) >= Lambda_v >= ell(v)+1,  so claim2 (N(s_k,v) >= n-1) holds whenever
  Lambda_v >= n-1  (guaranteed by ell(v) >= n-2, a fortiori by ell(v) >= n).

Also check the induction-fiber fact:  for v <=_L v' <=_L mu_v, is mu_{v'} = mu_v ?
(needed so the up-induction stays within a single dominant hull).
"""
import sys
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


def N(u, ds):
    acc = set()
    for w in ds:
        acc |= supp(u, tuple(w))
    return len(acc)


def run(n):
    perms = [tuple(p[i] for i in range(n)) for p in Permutation.all_permutations(n)]
    ds = {p: left_downset(Permutation(list(p)), n) for p in perms}
    Lam = {p: len(ds[p]) for p in perms}
    mu_of = {}
    for p in perms:
        vi = ~Permutation(list(p))
        mu = ~(vi.minimal_dominant_above())
        mu_of[p] = tuple(mu[i] for i in range(n))

    # (B+)  N(u,v) >= Lambda_v  for all u,v
    Bplus_ok = True
    worstBplus = (10.0, None)
    for v in perms:
        Lv = Lam[v]
        for u in perms:
            nv = N(u, ds[v])
            if nv < Lv:
                Bplus_ok = False
            ratio = nv / Lv
            if ratio < worstBplus[0]:
                worstBplus = (ratio, (u, v, nv, Lv))

    # fiber fact: v <=_L v' <=_L mu_v  =>  mu_{v'} = mu_v
    fiber_ok = True
    fiber_bad = None
    for v in perms:
        mu = mu_of[v]
        mu_ds = ds[mu]
        vset = ds[v]
        for vp in mu_ds:              # vp <=_L mu
            if vp in vset:            # and v <=_L vp  (i.e. vp in downset... careful: vset is DOWN from v)
                pass
        # we need vp with v <=_L vp <=_L mu: vp in [v,mu]_L.
        # v <=_L vp iff v in downset(vp) iff vset subset of ds[vp]
        for vp in mu_ds:
            if vset.issubset(ds[vp]):   # v <=_L vp
                if mu_of[vp] != mu:
                    fiber_ok = False
                    fiber_bad = (v, vp, mu, mu_of[vp])
                    break
        if not fiber_ok:
            break

    print(f"===== n={n} =====")
    r, info = worstBplus
    print(f"  (B+) N(u,v) >= Lambda_v : {Bplus_ok}   min ratio N/Lam = {r:.3f}  {info}")
    print(f"  (fiber) v<=_L v'<=_L mu_v => mu_v'=mu_v : {fiber_ok}"
          + (f"  BAD {fiber_bad}" if not fiber_ok else ""))


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
