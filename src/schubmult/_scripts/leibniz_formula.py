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


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Leibniz-style divided-difference sweep over RC-graphs")
    parser.add_argument("n", type=int, help="permutation size")
    parser.add_argument("k", type=int, help="max p to loop down from (1..k)")
    args = parser.parse_args(argv)

    n = args.n
    k = args.k
    start_time = time.time()

    from schubmult import Permutation, RCGraph, RCGraphRing, uncode

    rc_ring = RCGraphRing()
    tring = rc_ring @ rc_ring  # tensor-ring factory

    # single rc_ring accumulator for all Kogan-Kumar divdiff contributions
    kk_elem_total = rc_ring.zero

    # Accumulate tensor-ring element incrementally into `tensor_acc` using tring.zero

    perms = list(Permutation.all_permutations(n))

    for perm in perms:
        rc_iter = RCGraph.all_rc_graphs(perm, n - 1)

        for rc in rc_iter:
            for pval in range(k, 0, -1):  # start at k down to 1
                principal = RCGraph.principal_rc(uncode([pval]), n - 1)

                # tensor accumulator is local to this (rc, pval) check
                tensor_acc = tring.zero

                # termA: rc.divdiff_desc(1) x principal.dualpieri(...)
                try:
                    outs0 = rc.divdiff_desc(1)
                    outs1 = principal.dualpieri(Permutation([2, 1]), Permutation([2, 1]))
                    for a in outs0:
                        for b in outs1:
                            try:
                                tensor_acc = tensor_acc + tring((a, b[-1]))
                            except Exception:
                                traceback.print_exc()
                except Exception:
                    traceback.print_exc()

                # termB: rc (unchanged) x principal.divdiff_desc(1)
                try:
                    outs1b = principal.divdiff_desc(1)
                    for b in outs1b:
                        tensor_acc = tensor_acc + tring((rc, b))

                except Exception:
                    traceback.print_exc()

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
                # print(f"Accumulated {total_terms} pair-terms in {elapsed:.2f}s", file=sys.stderr)

                # Print the accumulated tensor-ring element (if we built one incrementally)
                try:
                    try:
                        sympy.pretty_print(tensor_acc)
                    except Exception:
                        traceback.print_exc()
                        print(repr(tensor_acc), file=sys.stderr)
                except Exception:
                    traceback.print_exc()
                    print("Failed to print tensor accumulator.", file=sys.stderr)

                # Print aggregated Kogan-Kumar element we built directly
                try:
                    print("\nAggregated Kogan-Kumar divdiff element (rc_ring):", file=sys.stderr)
                    if kk_elem_total is None:
                        print("Could not construct rc_ring.zero(); no aggregated element available.", file=sys.stderr)
                    else:
                        sympy.pretty_print(kk_elem_total)
                except Exception:
                    traceback.print_exc()
                    # fallback: show whatever diagnostic counts we have in acc
                    try:
                        print("Diagnostic accumulator:", file=sys.stderr)
                        print({repr(k): v for k, v in acc.items()})
                    except Exception:
                        traceback.print_exc()
                        print("No diagnostics available.", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
