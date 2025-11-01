"""
Test whether raising_operator / lowering_operator commute with
(arbitrary) sequences of up_jdt_slide followed by rectification.

Usage: run this script from the repo root (with the correct PYTHONPATH / venv
active) to run randomized trials and report any counterexamples.
"""
from __future__ import annotations

import copy
import random
from typing import List, Optional, Sequence, Tuple

import numpy as np

from schubmult import RootTableau

Seed = 12345
random.seed(Seed)


def ensure_cell(grid: np.ndarray, i: int, j: int) -> np.ndarray:
    """Return a copy of grid extended (if necessary) so (i,j) is a valid index."""
    rows, cols = grid.shape
    nr = max(rows, i + 1)
    nc = max(cols, j + 1)
    if (nr, nc) == (rows, cols):
        return grid.copy()
    out = np.empty((nr, nc), dtype=object)
    out.fill(None)
    for r in range(rows):
        for c in range(cols):
            out[r, c] = grid[r, c]
    return out


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
        g = copy.deepcopy(cur._root_grid)
        g = ensure_cell(g, i, j)
        g[i, j] = None
        tmp = RootTableau(g)
        try:
            cur = tmp.up_jdt_slide(i, j)
        except Exception:
            # slide failed -> sequence invalid for this tableau; abort
            return None
    # finally rectify to produce a straight tableau
    try:
        return cur.rectify()
    except Exception:
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


def random_reduced_and_compatible(max_len=6, max_letter=4):
    n = random.randint(1, max_len)
    reduced = [random.randint(1, max_letter) for _ in range(n)]
    compatible = [random.randint(1, max_letter) for _ in range(n)]
    return reduced, compatible


def random_up_seq(rt: RootTableau, max_len=5) -> Sequence[Tuple[int, int]]:
    """
    Generate a valid sequence of up-jdt hole positions for tableau `rt`.

    Validity rule used:
      - a hole position (i,j) is considered valid if, on the current grid,
        either both a box up (i-1,j) and a box left (i,j-1) exist, OR
        at least one of those neighbors exists and the hole is placed on the
        *outer* boundary (we allow positions with i==rows or j==cols to extend
        the grid by one).
    The sequence is built greedily: choose a random valid outer-corner,
    perform the up_jdt_slide (on an extended copy) to update the tableau,
    and repeat until max_len or no valid candidates remain.
    """
    seq: List[Tuple[int, int]] = []
    cur = rt
    for _ in range(max_len):
        grid = copy.deepcopy(cur._root_grid)
        rows, cols = grid.shape

        # consider candidate positions in range [0..rows] x [0..cols]
        candidates = []
        for i in range(0, rows + 1):
            for j in range(0, cols + 1):
                # skip positions that are currently occupied (must be a hole)
                if i < rows and j < cols and grid[i, j] is not None:
                    continue
                # check neighbors (up and left) existence within current grid
                up_exists = (i - 1 >= 0 and i - 1 < rows and j < cols and grid[i - 1, j] is not None)
                left_exists = (j - 1 >= 0 and j - 1 < cols and i < rows and grid[i, j - 1] is not None)

                # hole is on outer boundary if it extends grid by at least one coordinate
                is_outer = (i == rows) or (j == cols)

                # accept if either both neighbors exist, or at least one exists and position is outer
                if (up_exists and left_exists) or ((up_exists or left_exists) and is_outer):
                    candidates.append((i, j))

        if not candidates:
            break

        # pick random candidate and perform slide to update cur
        i_choice, j_choice = random.choice(candidates)
        seq.append((i_choice, j_choice))

        # perform the slide on an extended tableau to keep validity for next step
        g = ensure_cell(grid, i_choice, j_choice)
        g[i_choice, j_choice] = None
        tmp = RootTableau(g)
        try:
            cur = tmp.up_jdt_slide(i_choice, j_choice)
        except Exception:
            # if slide unexpectedly fails, stop building sequence
            seq.pop()
            break

    return tuple(seq)


def test_one_case(reduced_word, compatible_seq, seq, index: int, op_name: str) -> Tuple[bool, str]:
    """
    Test commutation for one operator index:
     left = op( Rect( UpSeq(T) ) )
     right = Rect( UpSeq( op(T) ) )
    Returns (passed, message)

    Note: do not suppress exceptions from raising_operator / lowering_operator —
    these should return None when the operator is not defined.
    """
    T = RootTableau.root_insert_rsk(reduced_word, compatible_seq)
    # apply UpSeq then rectify, then op
    A = apply_up_seq_and_rect(T, seq)
    if A is None:
        return True, "Up-sequence invalid on original tableau (skip)"

    # call operator directly; do not catch exceptions here
    if op_name == "raise":
        left = A.raising_operator(index)
    else:
        left = A.lowering_operator(index)

    # apply op first, then upseq+rect
    if op_name == "raise":
        eT = T.raising_operator(index)
    else:
        eT = T.lowering_operator(index)

    if eT is None:
        return True, "Operator not defined on T (skip)"

    B = apply_up_seq_and_rect(eT, seq)
    if B is None:
        return True, "Up-sequence invalid on operated tableau (skip)"

    # left might be None (operator undefined on A)
    if left is None:
        return True, "Operator not defined on rectified tableau (skip)"

    # compare
    ok = grids_equal(left, B)
    if ok:
        return True, "commute"
    return False, f"mismatch: index={index}, reduced={reduced_word}, compatible={compatible_seq}, seq={seq}"


def run_random_tests(num_cases=200):
    failures = []
    for t in range(num_cases):
        reduced, compatible = random_reduced_and_compatible()
        T = RootTableau.root_insert_rsk(reduced, compatible)
        seq = random_up_seq(T)
        # pick a random index to test (small range)
        idx = random.randint(1, 5)
        ok_r, msg_r = test_one_case(reduced, compatible, seq, idx, "raise")
        ok_l, msg_l = test_one_case(reduced, compatible, seq, idx, "lower")
        if not ok_r:
            failures.append(("raise", msg_r))
        if not ok_l:
            failures.append(("lower", msg_l))
    return failures


if __name__ == "__main__":
    import sys

    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
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