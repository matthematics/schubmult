"""Check N^up(u,v) <= Lambda(mdom v) * N(u,v)  and  C^{-1} <= Lambda(mdom v)."""
import sys
from schubmult import DSx, Permutation, uncode
from schubmult.symbolic import expand


def main(n):
    P = list(Permutation.all_permutations(n))
    inv = {w: frozenset(w.inversion_set) for w in P}
    def mdom(v): return Permutation(list(~uncode((~v).theta())))
    Lam = {}
    for v in P:
        md = mdom(v)
        if md not in Lam:
            Lam[md] = sum(1 for x in P if inv[x] <= inv[md])
    def N(u, v):
        pr = DSx(u) * DSx(v, "z")
        return len({k for k, c in pr.items() if expand(c) != 0})

    bad1 = bad2 = 0
    w1 = w2 = None
    for u in P:
        for v in P:
            md = mdom(v)
            Nv = N(u, v)
            Nup = N(u, md)
            L = Lam[md]
            if Nup > L * Nv:
                bad1 += 1
            r1 = Nup / (L * Nv)
            if w1 is None or r1 > w1[0]:
                w1 = (r1, u, v, Nup, L, Nv)
            # C^{-1} <= Lambda(mdom v)
            dd = min([x for x in P if mdom(x) == md], key=lambda w: len(inv[w]))
            Ndd = N(u, dd)
            if Nup > L * Ndd:
                bad2 += 1
            r2 = Nup / (L * Ndd)
            if w2 is None or r2 > w2[0]:
                w2 = (r2, u, v, Nup, L, Ndd)
    r, u, v, Nup, L, Nv = w1
    print(f"S_{n}: N^up <= Lambda(mdom v)*N(u,v): {bad1} failures; tightest ratio {r:.3f} u={u} v={v} Nup={Nup} L={L} N={Nv}")
    r, u, v, Nup, L, Ndd = w2
    print(f"     C^-1 <= Lambda(mdom v):          {bad2} failures; tightest ratio {r:.3f} u={u} v={v} Nup={Nup} L={L} Nddom={Ndd}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
