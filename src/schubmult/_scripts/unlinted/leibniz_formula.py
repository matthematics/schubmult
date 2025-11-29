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


def _kk_dd_to_rc_ring_elem(kk_dd, rc_ring):
    """
    Normalize kk_dd into an RCGraphRing element.
    Accepts:
      - an RCGraphRing element (returned as-is),
      - a mapping {rc_graph: coeff},
      - an iterable of rc_graphs (treated as coefficient 1 each).
    Returns an rc_ring element (or raises if conversion fails).
    """
    # already an rc_ring element?
    try:
        # crude duck-typing: try adding to zero to test
        zero = rc_ring.zero
        if hasattr(kk_dd, "__add__") and (kk_dd + zero) is not None:
            return kk_dd
    except Exception:
        # not an rc_ring element
        pass

    # mapping -> sum(coeff * rc_ring(key))
    if isinstance(kk_dd, dict):
        acc = rc_ring.zero
        for key, coeff in kk_dd.items():
            acc = acc + (coeff * rc_ring(key))
        return acc

    # iterable -> sum of rc_ring(items)
    acc = rc_ring.zero
    for item in kk_dd:
        acc = acc + (1 * rc_ring(item))
    return acc

def divdiff_pair(tring, rc1, rc2, n):
    from schubmult import Permutation, RCGraph, CrystalGraphTensor
    from schubmult import x, y
    s = Permutation([2,1])
    res = tring.zero
    
    dual = rc2.dualpieri(s, s)
    apiece = rc1.divdiff_desc(1)
    posdeg = False
    for a in apiece:
        for vlist, perm_list, b in dual:
            toadd = sympy.S.One
            for i in range(len(vlist)):
                for j in range(len(vlist[i])):
                    toadd *= x[i + 1] - y[vlist[i][j]]
                    posdeg = True
                    if vlist[i][j] != 1:
                        toadd = 0
                        
            b = b.weight_reflection(1)
            res = res + toadd * tring((a.resize(n), b.resize(n)))
            tt = CrystalGraphTensor(RCGraph.principal_rc(Permutation.w0(n), n).resize(n), b.resize(n))
            tt = tt.raising_operator(1)
            while tt is not None:
                res = res + toadd * tring((a.resize(n), tt.factors[1].resize(n)))
                tt = tt.raising_operator(1)
    if not posdeg:
        try:
            factor2, row = rc2.exchange_property(1, return_row=True)
            if row == 1:
                res += tring((rc1.resize(n), factor2.resize(n)))
        except ValueError:
            print("Caosubalsk")
            pass
    return res

def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Leibniz-style divided-difference sweep over RC-graphs")
    parser.add_argument("n", type=int, help="permutation size")
    parser.add_argument("k", type=int, help="max p to loop down from (1..k)")
    args = parser.parse_args(argv)

    n = args.n
    k = args.k
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

    perms = list(Permutation.all_permutations(n))

    for perm in perms:
        rc_iter = RCGraph.all_rc_graphs(perm, n)

        for rc in rc_iter:
            for pval in range(k, 0, -1):  # start at k down to 1

                principal = RCGraph.principal_rc(uncode([pval]), n)
                print("Expanding")
                pretty_print(tring((principal, rc)))
                # tensor accumulator is local to this (rc, pval) check
                tensor_acc = divdiff_pair(tring, principal, rc, n)

                # termA: rc.divdiff_desc(1) x principal.dualpieri(...)
                kk_elem_total = rc_ring.zero
                # s = Permutation([2,1])
                
                # try:
                #     pass
                    # outs1a0 = {rc}
                    # dualp = rc.dualpieri(s, s)
                    # #rc0, row = rc.exchange_property(1,return_row=True) 
                    # #outs1a = rc.exchange_property(1)
                    # if True:#row == 1:
                    #     outs1a = rc.divdiff_desc(1)
                    #     #dppi = principal.dualpieri(s,s)
                    #     principal_minus = RCGraph.principal_rc(uncode([0, pval - 1]), n)
                    #     outs1b = {principal_minus}#principal.divdiff_desc(1)
                    #     tensor_acc += sum([tring((a.resize(n), b.resize(n))) for a in outs1a0 for b in outs1b], start=tring.zero)
                    #     # print(outs0)
                    #     # print(outs1)
                    #     #tensor_acc = sum([tring((a.resize(n), b.resize(n))) for a in outs1a0 for b in outs1b], start=tring.zero)

                        
                    #     #for a in dual:
                    #     for vlist, perm_list, b in dualp:
                    #         toadd = S.One
                    #         for i in range(len(vlist)):
                    #             for j in range(len(vlist[i])):
                    #                 toadd *= x[i + 1] - y[vlist[i][j]]
                    #         #for a in outs1a:
                    #         #b = principal
                    #         b = b.weight_reflection(1)
                    #         tensor_acc = tensor_acc + toadd * tring((a.resize(n), b.resize(n)))
                    #         # w0 = RCGraph.principal_rc(Permutation.w0(n), n)
                    #         tt = CrystalGraphTensor(b.resize(n), principal_minus.resize(n))
                    #         tt = tt.raising_operator(1)
                    #         bb = b.resize(n).raising_operator(1)
                            
                    #         # if bb is not None:
                    #         #     tt = CrystalGraphTensor((a.resize(n),bb))
                    #         # while tt is not None:
                    #         #     tensor_acc = tensor_acc + toadd * tring((a.resize(n), tt.factors[1].resize(n)))
                    #         #     tt = tt.raising_operator(1)
                    #         while bb is not None:
                    #             tensor_acc = tensor_acc + toadd * tring((a.resize(n), bb.resize(n)))
                    #             bb = bb.raising_operator(1)
                    #     test1 = False
                    #     test2 = False
                    #     for (aa, bb) in tensor_acc:
                    #         if aa == principal_minus:
                    #             test1 = True
                    #         if bb == principal_minus:
                    #             test2 = True
                        #if not test1 and not test2:
                        
                        # if not test2 and test1:
                        #     tensor_acc += sum([tring((a.weight_reflection(1).resize(n), b.resize(n))) for a in outs1b0 for b in outs1a], start=tring.zero)
                        #     tensor_acc += sum([tring((b.weight_reflection(1).resize(n), a.resize(n))) for a in outs1a0 for b in outs1b], start=tring.zero)

                # except Exception:
                #     traceback.print_exc()

                # termB: rc (unchanged) x principal.divdiff_desc(1)
                # try:
                #     # outs1b = principal.divdiff_desc(1)
                #     outs1b = {principal}
                #     # dprc = rc.dualpieri(s,s)
                #     dprc = rc.divdiff_desc(1)#{((), (), rc)}
                #     #print(f"{dprc=}")
                #     #outs1a = {dprc}#{RCGraph([(),*a[-1].shiftup(1).normalize()]) for a in dprc}
                    
                #     for a in dprc:
                        
                #         for b in outs1b:
                #             #for a in outs1a:
                #             tensor_acc = tensor_acc + toadd * tring((a.weight_reflection(1).resize(n), b.resize(n)))

                # except Exception:
                #     traceback.print_exc()

                # Also compute the Kogan-Kumar insert/divdiff element for this rc and pval,
                # aggregate its divdiff_desc(1) into kk_acc so we can present it later.
                try:
                    # use pval (not global k) as the insertion parameter
                    kk_insert = rc.kogan_kumar_insert(pval, [1] * pval)
                    kk_dd = kk_insert.divdiff_desc(1)
                    
                    try:
                        kk_piece = _kk_dd_to_rc_ring_elem(kk_dd, rc_ring)
                        kk_elem_total = kk_elem_total + kk_piece
                    except Exception:
                        traceback.print_exc()

                except Exception:
                    traceback.print_exc()

                # elapsed = time.time() - start_time
                # print(f"Accumulated {total_terms} pair-terms in {elapsed:.2f}s", )

                # Print the accumulated tensor-ring element (if we built one incrementally)
                try:
                    try:
                        print("Wondrous tensor")
                        sympy.pretty_print(tensor_acc)
                    except Exception:
                        traceback.print_exc()
                        print(repr(tensor_acc), )
                except Exception:
                    traceback.print_exc()
                    print("Failed to print tensor accumulator.", )

                # Print aggregated Kogan-Kumar element we built directly
                try:
                    print("\nAggregated Kogan-Kumar divdiff element (rc_ring):", )
                    if kk_elem_total is None:
                        print("Could not construct rc_ring.zero(); no aggregated element available.", )
                    else:
                        sympy.pretty_print(kk_elem_total)
                except Exception:
                    traceback.print_exc()
                    # fallback: show whatever diagnostic counts we have in acc
                    try:
                        print("Diagnostic accumulator:", )
                        print({repr(k): v for k, v in acc.items()})
                    except Exception:
                        traceback.print_exc()
                        print("No diagnostics available.", )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
