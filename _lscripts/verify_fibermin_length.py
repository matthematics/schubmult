"""Check properties of the UNIQUE fiber minimum v0 = pi_down(mu):
  - length ell(v0) as mu ranges over dominant perms
  - is ell(v0) >= n  (except trivial small mu)?
  - does N(s_k, v0) >= n-1 for all fiber minima (except trivial)?
  - which minima violate, and are they exactly the "trivial" ones (mu with L small)?
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


def length(a, n):
    return sum(1 for i in range(n) for j in range(i + 1, n) if a[i] > a[j])


def run(n):
    perms = [tuple(p[i] for i in range(n)) for p in Permutation.all_permutations(n)]
    ds = {p: left_downset(Permutation(list(p)), n) for p in perms}
    mu_of = {}
    for p in perms:
        vi = ~Permutation(list(p))
        mu = ~(vi.minimal_dominant_above())
        mu_of[p] = tuple(mu[i] for i in range(n))
    from collections import defaultdict
    fibers = defaultdict(list)
    for p in perms:
        fibers[mu_of[p]].append(p)

    simple = []
    for k in range(n - 1):
        a = list(range(1, n + 1))
        a[k], a[k + 1] = a[k + 1], a[k]
        simple.append(tuple(a))

    rows = []
    for mu, elts in fibers.items():
        mins = [v for v in elts if not any((w != v and tuple(w) in ds[tuple(v)]) for w in elts)]
        v0 = mins[0]
        lv0 = length(v0, n)
        L = len([c for c in (~Permutation(list(mu))).code if c != 0])
        minN = min(N(s, ds[tuple(v0)]) for s in simple)
        rows.append((lv0, L, tuple(v0), tuple(mu), minN))

    rows.sort()
    print(f"===== n={n} =====")
    # how many fiber minima have ell(v0) < n, and do any NONtrivial (L>=2) violate N(s,v0)>=n-1?
    viol = [r for r in rows if r[4] < n - 1]
    short = [r for r in rows if r[0] < n]
    print(f"  #fibers={len(rows)}  #minima with ell(v0)<n: {len(short)}")
    print(f"  #minima with min_s N(s,v0) < n-1: {len(viol)}")
    print(f"  violators (ell(v0),L,v0,mu,minN):")
    for r in viol:
        print(f"     {r}")
    # relationship: is ell(v0) always >= something like n-1 for L>=2?
    print(f"  minima with L>=2 and ell(v0)<n: {[r for r in rows if r[1]>=2 and r[0]<n]}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
