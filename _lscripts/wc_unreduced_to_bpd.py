"""Prototype (user's tack):

For a WCGraph whose Hecke *column-insertion* tableau ``P`` gives an **unreduced**
word (``len(P.row_word) != perm.inv``):

1. Compute ``(P, Q) = WCGraph.hecke_invariant`` where
   ``Q`` is the set-valued *recording* tableau.
2. **Snap** ``Q`` to its minimum -> an ordinary (single-entry) recording tableau
   ``Qmin`` of the *same shape* as ``P``.
3. **Uninsert**: ``hecke_column_uninsert_rsk(P, Qmin)`` recovers the two-line
   array ``(recording_rows, insertion_word)``.  The insertion word here is the
   *real unreduced word* of the increasing tableau ``P`` (length = #boxes of P).
4. Feed the ``(letter, row)`` pairs to ``BPD.inverse_pop_op`` to try to build a
   BPD directly from that unreduced word.

Run:  conda activate schubmult_312 && python -O _lscripts/wc_unreduced_to_bpd.py
"""

from __future__ import annotations

import sys

import numpy as np

from schubmult import Permutation, uncode
from schubmult.combinatorics.bpd import BPD, TileType
from schubmult.combinatorics.increasing_tableau import IncreasingTableau
from schubmult.combinatorics.set_valued_tableau import SetValuedTableau
from schubmult.combinatorics.wc_graph import WCGraph


def snap_recording_to_min(Q):
    """Flatten a set-valued recording tableau to its minimum single entry per box."""
    cells = Q.cells if isinstance(Q, SetValuedTableau) else dict(Q)
    return {rc: (min(labels),) for rc, labels in cells.items()}


def find_unreduced_wcgraphs(w):
    n = len(w)
    out = []
    for wc in WCGraph.all_wc_graphs(w, n):
        P, Q = wc.hecke_invariant
        if len(P.row_word) != wc.perm.inv:
            out.append((wc, P, Q))
    return out


def try_uninsert_and_bpd(w, limit=6):
    n = len(w)
    ell = w.inv
    cands = find_unreduced_wcgraphs(w)
    print(f"\n=== w={w.trimcode} ell={ell}: {len(cands)} WCGraphs with UNREDUCED P ===")
    for wc, P, Q in cands[:limit]:
        pboxes = len(P.row_word)
        print(f"\n  WCGraph rows={tuple(tuple(r) for r in wc)}")
        print(f"    perm_word={wc.perm_word} compat={wc.compatible_sequence}")
        print(f"    P.row_word={P.row_word} (#boxes={pboxes}, ell={ell})")

        Qmin = snap_recording_to_min(Q)
        print(f"    Q (set-valued) shape={Q.shape if hasattr(Q,'shape') else '?'}  Qmin={Qmin}")

        try:
            recording, insertion = IncreasingTableau.hecke_column_uninsert_rsk(P, Qmin)
        except Exception as e:  # noqa: BLE001
            print(f"    uninsert FAILED: {type(e).__name__}: {e}")
            continue
        print(f"    uninsert -> recording(rows)={recording}  insertion(word)={insertion}")

        # verify the recovered word is a Hecke/ordinary word for w
        try:
            prod = Permutation.ref_product(*insertion)
            print(f"    ref_product(word)={prod.trimcode} inv={prod.inv} (==w? {prod == w}) len={len(insertion)}")
        except Exception as e:  # noqa: BLE001
            print(f"    ref_product FAILED: {type(e).__name__}: {e}")

        # build (letter,row) pairs for inverse_pop_op
        pairs = list(zip([int(x) for x in insertion], [int(r) for r in recording]))
        print(f"    (letter,row) pairs = {pairs}")

        base = BPD.rothe_bpd(Permutation([]), max(n, len(w)))
        try:
            bpd = base.inverse_pop_op(*pairs)
            blanks = int((bpd._grid == TileType.BLANK).sum())
            cross = int((bpd._grid == TileType.CROSS).sum())
            print(f"    inverse_pop_op -> perm={bpd.perm.trimcode} inv={bpd.perm.inv} "
                  f"blanks={blanks} cross={cross} valid={bpd.is_valid} reduced={bpd.is_reduced}")
            print(f"      perm==w? {bpd.perm == w}")
        except Exception as e:  # noqa: BLE001
            print(f"    inverse_pop_op FAILED: {type(e).__name__}: {e}")


if __name__ == "__main__":
    cases = [[1, 0, 2, 1], [2, 0, 2], [1, 0, 3]]
    if len(sys.argv) > 1:
        cases = [[int(c) for c in arg.split(",")] for arg in sys.argv[1:]]
    for c in cases:
        try_uninsert_and_bpd(uncode(c))
