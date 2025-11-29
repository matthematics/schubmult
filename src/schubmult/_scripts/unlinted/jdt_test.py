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



def apply_up_seq_and_rect(rt: RootTableau, seq: Sequence[Tuple[int, int]], index=None, lower=False) -> Optional[RootTableau]:
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
    if index is not None:
        if lower:
            cur = cur.lowering_operator(index)
        else:
            cur = cur.raising_operator(index)
    if cur is not None:
        return cur.rectify(randomized=True)
    return None


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
    if op_name == "raise":
        T2 = T.raising_operator(index)
        if T2 is None:
            ok = T.rc_graph.raising_operator(index) is None
            if ok:
                msg = "Raising annihilates both"
            else:
                msg = f"Raising annihilates left only, {T=} {T2=}"
        else:
            ok = (T2.rc_graph == T.rc_graph.raising_operator(index))
            if ok:
                msg = "Raising commutes"
            else:
                msg = f"Raising mismatch, {T=} {T.rc_graph=} vs {T2=} {T2.rc_graph=} {index=}"
    if op_name == "raiserectify":
        w2 = T.weight_tableau.raising_operator(index)

        seq = random_up_seq(T)
        if len(seq) == 0:
            return True, "empty up-seq, skipped"
        T2 = T.raising_operator(index)
        if T2 is not None:
            B = apply_up_seq_and_rect(T, seq, index=index)
            ok = (B == T2 and B.rc_graph.is_valid and B.perm == T.perm)
            if ok:
                msg = "Raising commutes"
            else:
                msg = f"Raising mismatch, {B.weight_tableau=} vs {w2=}, {B=} vs {T2=} {index=} {B.rc_graph=} {T2.rc_graph=} {B.perm=} vs {T.perm=} This is OK just simultaneous RC pairing {T.rc_graph=}"
        else:
            ok, msg = True, "Anniliate"
    if op_name == "lowerrectify":
        seq = random_up_seq(T)
        if len(seq) == 0:
            return True, "empty up-seq, skipped"
        T2 = T.lowering_operator(index)
        if T2 is not None:
            B = apply_up_seq_and_rect(T, seq, index=index, lower=True)
        else:
            B = None
        if T2 is not None:
            ok = (B == T2 and B.rc_graph.is_valid and B.perm == T.perm)
            if ok:
                msg = "Lowering commutes"
            else:
                msg = f"Lowering mismatch, {B.weight_tableau=} vs {w2=}, {B=} vs {T2=} {index=} {B.rc_graph=} {T2.rc_graph=} {B.perm=} vs {T.perm=} This is OK just simultaneous RC pairing {T.rc_graph=}"
        else:
            ok = B is None
            if ok:
                msg = "Lowering annihilates both"
            else:
                msg = f"Lowering annihilates left only, {T=} {B=}"
    elif op_name == "rectify":
        seq = random_up_seq(T)
        if len(seq) == 0:
            return True, "empty up-seq, skipped"
        B = apply_up_seq_and_rect(T, seq)

        # compare
        ok = T == B
        msg = "unique rectification"
    if ok:
        return True, msg
    # include RC/perm info if available
    # extra = f" rc={T} {B}"
    return False, msg




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

        outer_corners = tuple(cur.iter_outer_corners)
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
        rc = None
        for rc in rc_iter:
            break
        if rc is None:
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
        for idx in range(1, T.crystal_length()):
            ok_r, msg_r = test_one_case(T, idx, "raiserectify")
            if ok_r:
                logger.info("Case %d RECTIFY: OK — %s; idx=%s", t, msg_r, idx)
                pretty_print(T)
            else:
                logger.error("Case %d RECTIFY: FAIL — %s; idx=%s", t, msg_r, idx)
                pretty_print(T)
                # exit immediately on failure showing the failing case
                sys.exit(2)
    # if we reach here all cases passed or were skipped
    logger.info("All %d cases completed (no failing case encountered).", num_cases)
    return []

def run_complete_tests():
    """
    Run randomized tests and log the outcome of each trial (both raise and lower).
    Ensure the up-jdt sequence is actually valid for the sampled tableau: regenerate
    until apply_up_seq_and_rect(T, seq) succeeds (or give up after attempts).
    """
    # ensure logging is configured even if the module was imported
    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    t = 0
    for new_perm in _perms:
        if new_perm.inv == 0 or len(new_perm.trimcode) < 6:
            continue
        for rc in RCGraph.all_rc_graphs(new_perm, len(new_perm.trimcode)):
            t += 1
            T = RootTableau.from_rc_graph(rc)

            # build a valid up-sequence (retry until apply_up_seq_and_rect succeeds)
            # pick a random index to test (small range)
            idx = 0 #random.randint(1, 5)

            tests = ["raiserectify"]
            indexes = list(range(1, rc.crystal_length()))
            for test_op in tests:
                for idx in indexes:
                    ok_r, msg_r = test_one_case(T, idx, test_op)
                    if ok_r:
                        logger.info("Case %d %s: OK — %s; idx=%s", t, test_op, msg_r, idx)
                        pretty_print(T)
                    else:
                        logger.error("Case %d %s: FAIL — %s; idx=%s", t, test_op, msg_r, idx)
                        pretty_print(T)
                        # exit immediately on failure showing the failing case
                        sys.exit(2)
    # if we reach here all cases passed or were skipped
    logger.info("All %d cases completed (no failing case encountered).", t)
    return []



if __name__ == "__main__":
    import sys

    # simple logging setup so each test case is reported to stdout
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    perm_length = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    if perm_length > 15:
        raise ValueError("perm_length too large; may exhaust memory/time")
    N = int(sys.argv[2]) if len(sys.argv) > 2 else None

    _perms = Permutation.all_permutations(perm_length)
    print(f"Running {N} randomized tests (seed={Seed}) ...")
    if N is not None:
        fails = run_random_tests(N)
    else:
        fails = run_complete_tests()
    if not fails:
        print("No failures found (all tested cases either commute or were skipped).")
    else:
        print(f"Found {len(fails)} failures. Showing up to 10:")
        for k, f in enumerate(fails[:10], 1):
            print(k, f)
        # exit non-zero to indicate counterexamples found
        sys.exit(2)
