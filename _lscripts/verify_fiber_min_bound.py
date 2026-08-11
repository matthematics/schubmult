"""With the fiber = interval [v0,mu]_L known (Cambrian), find the SHARP bound at
the unique fiber minimum v0.  Candidates (only at v = v0 = pi_down(mu)):

  (C1)  N(u,mu) <= |[v0,mu]_L| * N(u,v0)            (fiber/class size)
  (C2)  N(u,mu) <= (Lambda(mu)/Lambda(v0)) * N(u,v0)  (P1p restricted to minima)
  (C3)  raw ratio  N(u,mu)/N(u,v0)  -- max, and compare to |[v0,mu]| and Lambda(mu)
Also: does |[v0,mu]_L| = prod(theta_i+1)/something?  print (Lambda(mu),Lambda(v0),
       |[v0,mu]|, L) per fiber.
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

    from collections import defaultdict
    fibers = defaultdict(list)
    for p in perms:
        fibers[mu_of[p]].append(p)

    worst = {"C1": (0.0, None), "C2": (0.0, None), "C3": (0.0, None)}

    def upd(name, val, info):
        if val > worst[name][0]:
            worst[name] = (val, info)

    for mu, elts in fibers.items():
        # unique minimum
        mins = [v for v in elts if not any((w != v and tuple(w) in ds[tuple(v)]) for w in elts)]
        v0 = mins[0]
        fiber_size = len(elts)                 # |[v0,mu]_L|
        Lmu = Lam[mu]
        Lv0 = Lam[tuple(v0)]
        theta = sorted((c for c in (~Permutation(list(mu))).code), reverse=True)
        # note: code(mu^{-1}) = theta(v^-1); parts:
        thf = [c for c in (~Permutation(list(mu))).code if c != 0]
        for u in perms:
            Nv0 = N(u, ds[tuple(v0)])
            if Nv0 == 0:
                continue
            Nmu = N(u, ds[mu])
            info = (u, tuple(v0), mu, fiber_size, Lmu, Lv0, tuple(sorted(thf, reverse=True)), Nv0, Nmu)
            upd("C1", Nmu / (fiber_size * Nv0), info)
            upd("C2", Nmu / ((Lmu / Lv0) * Nv0), info)
            upd("C3", Nmu / Nv0, info)

    print(f"===== n={n} =====")
    for name, label in (("C1", "|[v0,mu]| * N(u,v0)"), ("C2", "(Lam_mu/Lam_v0) N(u,v0)"), ("C3", "N(u,v0)")):
        r, info = worst[name]
        flag = "OK  " if (name != "C3" and r <= 1.0 + 1e-9) else ("    " if name == "C3" else "FAIL")
        print(f"  [{flag}] max N_mu/({label}) = {r:.4f}   {info}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
