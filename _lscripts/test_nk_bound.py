"""Aggregate statistics for candidate bounds on N(u, mu_v) vs N(u, v).

For all u,v in S_n:  N(u,w) = |union_{w'<=_L w} supp(S_u S_{w'})|,
mu_v = dom(v^{-1})^{-1}, k = #descents(v).
Report the worst (max) value of several candidate ratios.
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


def run(n):
    perms = Permutation.all_permutations(n)
    stats = {
        "raw   N_mu/N_v": (0.0, None),
        "clean /(n^k Nv)": (0.0, None),
        "quad  /(n^k Nv^2)": (0.0, None),
        "quadplain /Nv^2": (0.0, None),
        "vs interval ratio": (0.0, None),
    }

    def upd(name, val, info):
        if val > stats[name][0]:
            stats[name] = (val, info)

    for v in perms:
        k = len(v.descents())
        mu = ~((~v).minimal_dominant_above())
        v_ds = left_weak_downset(v, n)
        mu_ds = left_weak_downset(mu, n)
        assert tuple(v[i] for i in range(n)) in mu_ds
        nk = n ** k
        intratio = len(mu_ds) / len(v_ds)
        for u in perms:
            uu = tuple(u[i] for i in range(n))
            Nv = N(uu, v_ds)
            Nmu = N(uu, mu_ds)
            if Nv == 0:
                continue
            info = (tuple(u), tuple(v), tuple(mu), k, Nv, Nmu, nk)
            upd("raw   N_mu/N_v", Nmu / Nv, info)
            upd("clean /(n^k Nv)", Nmu / (nk * Nv), info)
            upd("quad  /(n^k Nv^2)", Nmu / (nk * Nv * Nv), info)
            upd("quadplain /Nv^2", Nmu / (Nv * Nv), info)
            upd("vs interval ratio", (Nmu / Nv) / intratio, info)
    print(f"===== n={n} =====")
    for name, (val, info) in stats.items():
        print(f"  max {name:20s} = {val:.4f}   {info}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5"])):
        run(n)
