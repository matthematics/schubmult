"""[e, mdom v]_L is a union of whole Tamari classes, giving an exact decomposition
   Lambda(mdom v) = sum of |C| over classes C with top <=_L mdom v.
Also: how many dominants lie below mdom v, vs Lambda(ddom v)?"""

import sys

from schubmult import Permutation, uncode


def main(n):
    P = list(Permutation.all_permutations(n))
    inv = {w: frozenset(w.inversion_set) for w in P}

    def mdom(v):
        return Permutation(list(~uncode((~v).theta())))

    cls = {}
    for v in P:
        cls.setdefault(mdom(v), []).append(v)

    badunion = 0
    badexact = 0
    baddom = 0
    for md, cl in cls.items():
        below = [x for x in P if inv[x] <= inv[md]]
        # union of whole classes?
        tops = {mdom(x) for x in below}
        for t in tops:
            if not all(inv[y] <= inv[md] for y in cls[t]):
                badunion += 1
        if len(below) != sum(len(cls[t]) for t in tops):
            badexact += 1
        dd = min(cl, key=lambda w: len(inv[w]))
        lam_dd = sum(1 for x in P if inv[x] <= inv[dd])
        if len(tops) > lam_dd:
            baddom += 1
            if baddom <= 4:
                print(f"   #classes below {md} = {len(tops)} > Lambda(ddom={dd}) = {lam_dd}")
    print(f"S_{n}: {len(cls)} classes")
    print(f"  [e,mdom]_L is a union of whole classes : {badunion} failures")
    print(f"  exact decomposition Lambda = sum |C|   : {badexact} failures")
    print(f"  #classes below mdom <= Lambda(ddom)    : {baddom} failures")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
