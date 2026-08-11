"""Test the LINEAR bound  N(u,mu_v) <= n^param * N(u,v)  for param in {k, L}.
   k = #descents(v);  L = #nonzero parts of theta(v^{-1})  (the column count).
"""
import sys
from schubmult import Permutation
from schubmult.rings.schubert.schubert_ring import Sx


def left_weak_downset(perm, n):
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


def N(u, downset):
    acc = set()
    for w in downset:
        acc |= supp(u, w)
    return len(acc)


def run(n, u_sample=None):
    perms = Permutation.all_permutations(n)
    best = {"k": (0.0, None), "L": (0.0, None), "D": (0.0, None), "D2": (0.0, None)}
    mismatch_LD = []  # v where L != maxdesc(v)

    def upd(name, val, info):
        if val > best[name][0]:
            best[name] = (val, info)

    for v in perms:
        vi = ~v
        k = len(v.descents())
        L = len([c for c in vi.theta() if c != 0])
        # C = number of nonzero parts of code(v^{-1})  (should equal L)
        C = len([c for c in vi.code if c != 0])
        # max descent of v^{-1} = number of x-variables needed for S_{v^{-1}} (1-indexed)
        desci = vi.descents()
        D = (max(desci) + 1) if desci else 0
        if L != C:
            mismatch_LD.append((tuple(v), L, C))
        mu = ~(vi.minimal_dominant_above())
        v_ds = left_weak_downset(v, n)
        mu_ds = left_weak_downset(mu, n)
        us = perms if u_sample is None else u_sample
        for u in us:
            uu = tuple(u[i] for i in range(n))
            Nv = N(uu, v_ds)
            if Nv == 0:
                continue
            Nmu = N(uu, mu_ds)
            info = (tuple(u), tuple(v), tuple(mu), k, L, D, Nv, Nmu)
            upd("k", Nmu / (n ** k * Nv), info)
            upd("L", Nmu / (n ** L * Nv), info)
            upd("D", Nmu / (n ** D * Nv), info)
            upd("D2", Nmu / (n ** D * Nv * Nv), info)
    print(f"===== n={n} =====")
    print(f"  L != #nonzero code(v^-1) count: {len(mismatch_LD)}"
          + (f"  e.g. {mismatch_LD[:3]}" if mismatch_LD else "  (L == #nonzero code(v^-1) always)"))
    for name in ("k", "L", "D", "D2"):
        r, info = best[name]
        flag = "OK  " if r <= 1.0 + 1e-9 else "FAIL"
        label = {"k": "n^k  Nv", "L": "n^L  Nv", "D": "n^D  Nv", "D2": "n^D Nv^2"}[name]
        print(f"  [{flag}] max N_mu/({label}) = {r:.4f}   {info}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
