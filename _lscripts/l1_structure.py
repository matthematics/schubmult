#!/usr/bin/env python3
"""Structure of the L=1 case of the dominant hull lemma, and the order relation
between v and mdom v.  Permutation arithmetic only.
"""

import sys

from schubmult.combinatorics.permutation import Permutation


def L_of(v):
    return len([t for t in v.theta() if t != 0])


def main(n):
    print(f"n={n}")
    bad_le = 0
    ones = []
    for v in Permutation.all_permutations(n):
        vi = ~v
        b = v.mul_dominant()
        a = v.mul_sortable()
        if not v.weak_order_leq(b) if hasattr(v, "weak_order_leq") else False:
            bad_le += 1
        if L_of(vi) == 1:
            ones.append((v, a, b, vi.theta()))
    print(f"   v with L(theta(v^-1))=1 : {len(ones)}")
    for v, a, b, th in ones:
        print(f"      v={tuple(v)!s:<16} code(v)={tuple(v.code)!s:<14} theta(v^-1)={tuple(t for t in th if t)!s:<6} ddom={tuple(a)!s:<14} mdom={tuple(b)!s}")


if __name__ == "__main__":
    for n in range(3, (int(sys.argv[1]) if len(sys.argv) > 1 else 4) + 1):
        main(n)
        print()
