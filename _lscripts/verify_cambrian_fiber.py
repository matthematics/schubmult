"""Verify the Cambrian / lattice-congruence structure of the dominant hull
   mu_v := dom(v^{-1})^{-1}  =  pi^uparrow(v)   (Reading 2005).

Claims to check exhaustively:
  (I)  each fiber F_mu = {v : mu_v = mu} has a UNIQUE <=_L-minimum v0.
  (II) F_mu is exactly the interval [v0, mu]_L  (congruence class = interval / convex).
  (III) pi^uparrow order-preserving: v <=_L v'  =>  mu_v <=_L mu_v'.
  (IV) idempotent: mu_{mu_v} = mu_v; and mu_v is dominant (132-avoiding).
  (V) #distinct mu (= #dominant perms) = Catalan(n).
"""
import sys
from math import comb
from schubmult import Permutation


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


def leq_L(a, b, n, ds_cache):
    # a <=_L b  iff  a in downset(b)
    return tuple(a) in ds_cache[tuple(b)]


def is_dominant(a, n):
    # 132-avoiding
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if a[i] < a[k] < a[j]:
                    return False
    return True


def run(n):
    perms = [tuple(p[i] for i in range(n)) for p in Permutation.all_permutations(n)]
    ds = {p: left_downset(Permutation(list(p)), n) for p in perms}
    mu_of = {}
    for p in perms:
        vi = ~Permutation(list(p))
        mu = ~(vi.minimal_dominant_above())
        mu_of[p] = tuple(mu[i] for i in range(n))

    # fibers
    from collections import defaultdict
    fibers = defaultdict(list)
    for p in perms:
        fibers[mu_of[p]].append(p)

    unique_min = True
    interval_ok = True
    for mu, elts in fibers.items():
        elset = set(elts)
        # minimal elements
        mins = [v for v in elts if not any((w != v and tuple(w) in ds[tuple(v)]) for w in elts)]
        if len(mins) != 1:
            unique_min = False
            v0 = None
        else:
            v0 = mins[0]
            # interval check: [v0,mu]_L = {w : v0<=_L w <=_L mu} should equal elset
            interval = set()
            for w in perms:
                if tuple(v0) in ds[tuple(w)] and tuple(w) in ds[mu]:  # v0<=_L w<=_L mu
                    interval.add(w)
            if interval != elset:
                interval_ok = False

    # (III) order preserving
    order_ok = True
    for v in perms:
        for vp in ds[tuple(v)]:   # vp <=_L v
            if mu_of[tuple(vp)] not in ds[mu_of[tuple(v)]]:  # mu_vp <=_L mu_v ?
                order_ok = False
                break
        if not order_ok:
            break

    # (IV) idempotent + dominant
    idem_ok = all(mu_of[mu_of[p]] == mu_of[p] for p in perms)
    dom_ok = all(is_dominant(mu_of[p], n) for p in perms)

    # (V) catalan
    ndist = len(fibers)
    cat = comb(2 * n, n) // (n + 1)

    print(f"===== n={n} =====")
    print(f"  (I) unique fiber minimum: {unique_min}")
    print(f"  (II) fiber = interval [v0,mu]_L: {interval_ok}")
    print(f"  (III) pi^up order-preserving: {order_ok}")
    print(f"  (IV) idempotent: {idem_ok};  image dominant(132-avoid): {dom_ok}")
    print(f"  (V) #distinct mu = {ndist};  Catalan({n}) = {cat};  match: {ndist == cat}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["3", "4", "5"])):
        run(n)
