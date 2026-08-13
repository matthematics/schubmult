"""Loss factors for restating the complexity bound in Ntilde or N(u,v).

A1 = N(u,ddom v), Am = N(u,mdom v) = N^up, m = class size, Ntil = class average.
  (i)  N^up / Ntil      <= m           (trivial: Ntil >= Am/m)
  (ii) Ntil  / N(u,v)   <= (1+(m-1)/C)/m <= 1/C
  (iii) N^up / N(u,v)   <= 1/C
Measure how big these really get.
"""
import sys
from schubmult import DSx, Permutation, uncode
from schubmult.symbolic import expand


def main(n):
    P = list(Permutation.all_permutations(n))

    def mdom(v):
        return Permutation(list(~uncode((~v).theta())))

    inv = {w: frozenset(w.inversion_set) for w in P}
    cls = {}
    for v in P:
        cls.setdefault(mdom(v), []).append(v)
    bots = {}
    for md, cl in cls.items():
        mn = min(cl, key=lambda w: len(inv[w]))
        bots[md] = mn

    lam = {md: sum(1 for x in P if inv[x] <= inv[md]) for md in cls}

    def N(u, v):
        pr = DSx(u) * DSx(v, "z")
        return len({k for k, c in pr.items() if expand(c) != 0})

    w_up_til = w_til_N = w_up_N = None
    w_m = 0
    for u in P:
        for md, cl in cls.items():
            vals = {v: N(u, v) for v in cl}
            A1, Am = vals[bots[md]], vals[md]
            m = len(cl)
            til = sum(vals.values()) / m
            C = A1 / Am
            r1 = Am / til
            if w_up_til is None or r1 > w_up_til[0]:
                w_up_til = (r1, u, md, m, A1, Am, til)
            for v, Nv in vals.items():
                r2 = til / Nv
                r3 = Am / Nv
                if w_til_N is None or r2 > w_til_N[0]:
                    w_til_N = (r2, u, v, m, A1, Am, til, 1 / C)
                if w_up_N is None or r3 > w_up_N[0]:
                    w_up_N = (r3, u, v, m, A1, Am, til, 1 / C)
            w_m = max(w_m, m)

    print(f"S_{n}: {len(P)} perms, {len(cls)} classes, max class size m={w_m}, max Lambda(mdom)={max(lam.values())}")
    r, u, md, m, A1, Am, til = w_up_til
    print(f"  max N^up/Ntilde = {r:.4f}  (bound m={m})   u={u} class top={md} A1={A1} Am={Am} Ntil={til:.3f}")
    r, u, v, m, A1, Am, til, invC = w_til_N
    print(f"  max Ntilde/N(u,v) = {r:.4f}  (bound 1/C={invC:.3f})  u={u} v={v} m={m} A1={A1} Am={Am}")
    r, u, v, m, A1, Am, til, invC = w_up_N
    print(f"  max N^up/N(u,v)  = {r:.4f}  (bound 1/C={invC:.3f})  u={u} v={v} m={m} A1={A1} Am={Am}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
