"""Sufficient conditions to prove (P1)  N(u,mu_v) <= Lambda(mu_v) N(u,v)
via a union bound over w' <=_L mu_v:

  N(u,mu_v) = | U_{w'<=_L mu_v} supp(S_u S_{w'}) |.

Candidate per-term bounds (if any holds for ALL w' <=_L mu_v, union bound closes it):
  (Q1)  |supp(S_u S_{w'})|                 <= N(u,v)        for all w'<=_L mu_v
  (Q2)  |supp(S_u S_{w'})|                 <= N(u,v)        only for w' with l(w')<=l(v)? 
Also directly measure the union "efficiency": how many DISTINCT supp sets, and their max size.
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
        acc |= supp(u, w)
    return len(acc)


def run(n):
    perms = Permutation.all_permutations(n)
    worstQ1 = (0.0, None)     # max over (u,v,w') of |supp(u,w')| / N(u,v)
    for v in perms:
        vi = ~v
        mu = ~(vi.minimal_dominant_above())
        lv = v.inv  # length of v
        v_ds = left_downset(v, n)
        mu_ds = left_downset(mu, n)
        for u in perms:
            uu = tuple(u[i] for i in range(n))
            Nv = N(uu, v_ds)
            if Nv == 0:
                continue
            for wp in mu_ds:
                s = len(supp(uu, wp))
                r = s / Nv
                if r > worstQ1[0]:
                    worstQ1 = (r, (tuple(u), tuple(v), tuple(mu), wp, s, Nv))
    r, info = worstQ1
    flag = "OK  " if r <= 1.0 + 1e-9 else "FAIL"
    print(f"n={n}: (Q1) max_w' |supp(S_u S_w')| / N(u,v) = {r:.4f}  [{flag}]  {info}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
