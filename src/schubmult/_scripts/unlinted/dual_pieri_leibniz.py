"""
Leibniz-style sweep over RC-graphs producing a tensor in (RCGraphRing @ RCGraphRing).

Usage:
    python leibniz_formula.py N K

All exceptions print a full traceback to stderr.
"""

from __future__ import annotations

import argparse
import collections
import sys
import time
import traceback
import sympy
from typing import List



# def divdiff_pair(tring, rc1, rc2, n):
#     from schubmult import Permutation, RCGraph, CrystalGraphTensor
#     from schubmult import x, y
#     s = Permutation([2,1])
#     res = tring.zero
    
#     dual = rc2.dualpieri(s, s)
#     apiece = rc1.divdiff_desc(1)
#     posdeg = False
#     for a in apiece:
#         for vlist, perm_list, b in dual:
#             if 
#     factor2 = rc2.divdiff_desc(1)
#     for b in factor2:
#         res += tring((rc1.resize(n), b.resize(n)))
#     return res

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Leibniz-style divided-difference sweep over RC-graphs")
    parser.add_argument("n", type=int, help="permutation size")
    #parser.add_argument("k", type=int, help="max p to loop down from (1..k)")
    args = parser.parse_args(argv)

    n = args.n
    #k = args.k
    start_time = time.time()

    from schubmult import Permutation, RCGraph, RCGraphRing, uncode, CrystalGraphTensor
    from schubmult import x, y
    from symengine import S
    from sympy import pretty_print

    rc_ring = RCGraphRing()
    tring = rc_ring @ rc_ring  # tensor-ring factory

    # single rc_ring accumulator for all Kogan-Kumar divdiff contributions
    # kk_elem_total = rc_ring.zero

    # Accumulate tensor-ring element incrementally into `tensor_acc` using tring.zero
    s = Permutation([2,1])
    perms = list(Permutation.all_permutations(n))

    dominant_perms = {perm for perm in perms if perm.minimal_dominant_above() == perm}

    
    for perm in perms:

        print(f"{perm.trimcode=}")

        rc_iter = RCGraph.all_rc_graphs(perm, n)
        
        for rc in rc_iter:
            for dom in dominant_perms:
                print(f"{dom.trimcode=}")
                print("RC:")
                pretty_print(rc)
                # left dominant, right arbitrary, s1 div diff
                pret = rc.dualpieri(dom,dom)
                for dual in pret:
                    print("dualpieri:")
                    pretty_print((dual[0], dual[-1]))
        print("==============================================================")
                
    print("paint")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
