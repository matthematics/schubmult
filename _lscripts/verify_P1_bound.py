"""Test candidate clean bounds that would PROVE  N(u,mu_v) <= n^L N(u,v):

  Lambda(mu_v) := |[e,mu_v]_L|   (proven <= prod(theta_i+1) <= n^L).

  (P1)  N(u,mu_v) <= Lambda(mu_v) * N(u,v)
  (P1') N(u,mu_v) <= (Lambda(mu_v)/Lambda(v)) * N(u,v)     Lambda(v):=|[e,v]_L|
  (P2)  N(u,mu_v) <= n^L * N(u,v)                          (target)

If (P1) holds exhaustively, then (P2) follows from Lambda(mu_v) <= n^L.
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
    worst = {"P1": (0.0, None), "P1p": (0.0, None), "P2": (0.0, None)}

    def upd(name, val, info):
        if val > worst[name][0]:
            worst[name] = (val, info)

    for v in perms:
        vi = ~v
        theta = sorted((c for c in vi.code), reverse=True)
        L = len([c for c in theta if c != 0])
        mu = ~(vi.minimal_dominant_above())
        v_ds = left_downset(v, n)
        mu_ds = left_downset(mu, n)
        Lam_mu = len(mu_ds)          # |[e,mu]_L|
        Lam_v = len(v_ds)            # |[e,v]_L|
        for u in perms:
            uu = tuple(u[i] for i in range(n))
            Nv = N(uu, v_ds)
            if Nv == 0:
                continue
            Nmu = N(uu, mu_ds)
            info = (tuple(u), tuple(v), tuple(mu), L, Lam_v, Lam_mu, Nv, Nmu)
            upd("P1", Nmu / (Lam_mu * Nv), info)
            upd("P1p", Nmu / ((Lam_mu / Lam_v) * Nv), info)
            upd("P2", Nmu / (n ** L * Nv), info)
    print(f"===== n={n} =====")
    for name, label in (("P1", "Lam_mu * Nv"), ("P1p", "(Lam_mu/Lam_v) * Nv"), ("P2", "n^L * Nv")):
        r, info = worst[name]
        flag = "OK  " if r <= 1.0 + 1e-9 else "FAIL"
        print(f"  [{flag}] max N_mu / ({label}) = {r:.4f}   {info}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
