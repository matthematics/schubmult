"""
Test whether raising_operator / lowering_operator commute with
(arbitrary) sequences of up_jdt_slide followed by rectification.

Usage: run this script from the repo root (with the correct PYTHONPATH / venv
active) to run randomized trials and report any counterexamples.
"""
from __future__ import annotations

import copy
import logging
import random
import sys
from typing import Any, List, Optional, Sequence, Tuple

import numpy as np
from sympy import pretty_print

from schubmult import Permutation, RCGraph, RootTableau

Seed = 12345
random.seed(Seed)

logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logger.addHandler(handler)
logger.setLevel(logging.INFO)
logger.propagate = False



def apply_up_seq_and_rect(rt: RootTableau, seq: Sequence[Tuple[int, int]]) -> Optional[RootTableau]:
    """
    For a given RootTableau (assumed straight), apply a sequence of "create hole
    at (i,j) then up_jdt_slide(i,j)" operations and finally rectify().

    Each (i,j) may extend the tableau (a hole on an outer corner). We create
    that hole by extending the underlying grid and setting that cell to None,
    then calling up_jdt_slide on the new tableau.
    """
    cur = rt
    for (i, j) in seq:
        # create a copy of current grid extended to include (i,j) and set hole
        logger.debug("Applying up-jdt slide at (%d,%d):", i, j)
        # pretty_print(cur)
        cur = cur.up_jdt_slide(i, j)
    logger.debug("After applying up-seq:")
    # pretty_print(cur)
    return cur.rectify()


def grids_equal(a: RootTableau, b: RootTableau) -> bool:
    """Compare two RootTableau by their underlying root grids (shape + cell equality)."""
    A = a._root_grid
    B = b._root_grid
    if A.shape != B.shape:
        return False
    for i, j in np.ndindex(A.shape):
        if A[i, j] is None and B[i, j] is None:
            continue
        if (A[i, j] is None) != (B[i, j] is None):
            return False
        if A[i, j] != B[i, j]:
            return False
    return True

_perms = tuple()
used_rcs = set()

# def sample_tableau_from_perms() -> tuple[RootTableau, Any]:
#     """
#     Pick a random permutation from _perms, pick a random RCGraph from
#     RCGraph.all_rc_graphs(perm, len(perm.trimcode)), and construct a RootTableau
#     via RootTableau.from_rc_graph(rc). Returns (RootTableau, rc) for logging.
#     """
#     if not _perms:
#         raise RuntimeError("No permutations available in _perms")

#     perm = random.choice(_perms)
#     # generate all rc graphs for this perm and pick one at random
#     all_rc = tuple(rc for rc in RCGraph.all_rc_graphs(perm, len(perm.trimcode)) if rc not in used_rcs)
#     if not all_rc:
#         _perms.remove(perm)
#         return sample_tableau_from_perms()
#     rc = random.choice(all_rc)
#     used_rcs.add(rc)
#     T = RootTableau.from_rc_graph(rc)
#     return T, rc


def test_one_case(T: RootTableau, index: int, op_name: str, rc=None) -> Tuple[bool, str]:
    """
    Test commutation for one operator index:
     left = op( Rect( UpSeq(T) ) )
     right = Rect( UpSeq( op(T) ) )
    Returns (passed, message)

    Note: do not suppress exceptions from raising_operator / lowering_operator —
    these should return None when the operator is not defined.
    """
    # apply UpSeq then rectify, then op


    # call operator directly; do not catch exceptions here
    seq = random_up_seq(T)
    if len(seq) == 0:
        return True, "empty up-seq, skipped"
    B = apply_up_seq_and_rect(T, seq)

    # compare
    ok = T.rc_graph == B.rc_graph
    if ok:
        return True, f"unique rectification, {len(seq)=}"
    # include RC/perm info if available
    extra = f" rc={T} {B}"
    return False, f"mismatch: index={index}, seq={seq}{extra}"


    

def random_up_seq(rt: RootTableau, max_len=25) -> Sequence[Tuple[int, int]]:
    """
    Generate a valid sequence of up-jdt hole positions for tableau `rt`.

    Each chosen (i,j) satisfies the outer-corner predicate above w.r.t. the
    current tableau. After choosing a hole we perform the up_jdt_slide to
    update the tableau so subsequent choices remain valid.
    """
    seq: List[Tuple[int, int]] = []
    cur = rt
    seq = []
    for _ in range(max_len):
        
        outer_corners = tuple(cur.iter_outer_corners())
        if not outer_corners:
            break
        box  = random.choice(outer_corners)
        seq = seq + [box]
        cur = cur.up_jdt_slide(*box)

    logger.debug(f"Generated sequence {seq=}")
    return tuple(seq)


def run_random_tests(num_cases=200):
    """
    Run randomized tests and log the outcome of each trial (both raise and lower).
    Ensure the up-jdt sequence is actually valid for the sampled tableau: regenerate
    until apply_up_seq_and_rect(T, seq) succeeds (or give up after attempts).
    """
    # ensure logging is configured even if the module was imported
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    perm_iter = iter(_perms)
    new_perm = next(perm_iter)
    while new_perm.inv == 0:
        new_perm = next(perm_iter)
    rc_iter = iter(RCGraph.all_rc_graphs(new_perm, len(new_perm.trimcode)))
    for t in range(1, num_cases + 1):

        try:
            rc = next(rc_iter)
        except StopIteration:
            try:
                new_perm = next(perm_iter)
            except StopIteration:
                logger.info("Exhausted all permutations after %d cases.", t - 1)
                break
            rc_iter = iter(RCGraph.all_rc_graphs(new_perm, len(new_perm.trimcode)))
            continue
        T = RootTableau.from_rc_graph(rc)

        # build a valid up-sequence (retry until apply_up_seq_and_rect succeeds)
        # pick a random index to test (small range)
        idx = 0 #random.randint(1, 5)

        ok_r, msg_r = test_one_case(T, idx, "rectify", rc=rc)
        if ok_r:
            logger.info("Case %d RECTIFY: OK — %s; idx=%s rc=%s", t, msg_r, idx, rc)
        else:
            logger.error("Case %d RECTIFY: FAIL — %s; idx=%s rc=%s", t, msg_r, idx, rc)
            # exit immediately on failure showing the failing case
            sys.exit(2)
    # if we reach here all cases passed or were skipped
    logger.info("All %d cases completed (no failing case encountered).", num_cases)
    return []


if __name__ == "__main__":
    import sys

    # simple logging setup so each test case is reported to stdout
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    perm_length = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    _perms = Permutation.all_permutations(perm_length)
    print(f"Running {N} randomized tests (seed={Seed}) ...")
    fails = run_random_tests(N)
    if not fails:
        print("No failures found (all tested cases either commute or were skipped).")
    else:
        print(f"Found {len(fails)} failures. Showing up to 10:")
        for k, f in enumerate(fails[:10], 1):
            print(k, f)
        # exit non-zero to indicate counterexamples found
        sys.exit(2)