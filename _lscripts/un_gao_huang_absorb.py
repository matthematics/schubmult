"""Prototype: un-Gao-Huang an *unreduced* BPD into a WCGraph.

Idea (from the analysis of ``BPD.pop_op``):

* ``pop_op`` is *inversion-driven*: each call removes exactly one inversion of
  the permutation and returns ``(new_bpd, (reflection, row))``.  The loop in
  ``as_reduced_compatible`` runs ``while perm.inv > 0``, so it pops exactly
  ``ell(w)`` times.

* On a *reduced* BPD ``#blanks == ell(w)``, so inversions and blanks are
  consumed in lockstep and every blank is popped.

* On an *unreduced* BPD there are ``#blanks - ell(w)`` *excess* blanks.  Each
  excess blank is a *second* crossing of a pair of pipes -- a repeated letter,
  invisible to ``perm.inv``.  Two pathologies result:

    - the driving mark lands on an excess blank -> ``pop_op`` cycles (same grid
      returns);  this is a **K-Knuth absorption** (Hecke ``alpha = 0``).
    - the excess blank sits aside the genuine inversions -> the loop stops early
      at ``perm.inv == 0`` leaving surplus blanks un-popped.

This prototype replaces the *inversion-driven* loop with a *blank-driven* loop
that detects absorptions:  pop until **no blanks remain**, and whenever a pop
step fails to lower ``perm.inv`` (cycle / no-progress), record that step as an
**absorption** and neutralise the offending blank so the loop can proceed.

Run:  conda activate schubmult_312 && python -O _lscripts/un_gao_huang_absorb.py
"""

from __future__ import annotations

import sys

import numpy as np

from schubmult import Permutation, uncode
from schubmult.combinatorics.bpd import BPD, TileType
from schubmult.combinatorics.wc_graph import WCGraph


def _num_blanks(bpd: BPD) -> int:
    return int((bpd._grid == TileType.BLANK).sum())


def _grid_key(bpd: BPD):
    return bpd._grid.tobytes(), bpd._grid.shape


def _driving_blank(bpd: BPD):
    """Return (row, col) of the blank ``pop_op`` would drive: rightmost blank in
    the first blank row."""
    blanks = np.argwhere(bpd._grid == TileType.BLANK)
    if len(blanks) == 0:
        return None
    r = blanks[0, 0]
    y = np.max(blanks[blanks[:, 0] == r, 1])
    return int(r), int(y)


def _neutralise_blank(bpd: BPD, r: int, c: int) -> BPD | None:
    """Try every tile substitution at (r, c); return the first that stays valid
    and strictly reduces the blank count (with unchanged permutation)."""
    inv0 = bpd.perm.inv
    for t in (TileType.CROSS, TileType.VERT, TileType.HORIZ, TileType.ELBOW_NW, TileType.ELBOW_SE):
        new = bpd.copy()
        new._grid[r, c] = t
        new.rebuild()
        if new.is_valid and _num_blanks(new) < _num_blanks(bpd) and new.perm.inv == inv0:
            return new
    return None


def diagnose(bpd: BPD):
    """Diagnostic: what does pop_op return on this (possibly unreduced) BPD, and
    which tile substitutions at the driving blank keep validity?"""
    inv0 = bpd.perm.inv
    print(f"    diagnose: perm={bpd.perm.trimcode} inv={inv0} blanks={_num_blanks(bpd)} reduced={bpd.is_reduced}")
    db = _driving_blank(bpd)
    print(f"    driving blank = {db}")
    try:
        popped, (refl, row) = bpd.pop_op()
        print(f"    pop_op -> letter=({refl},{row}) new_inv={popped.perm.inv} (delta {popped.perm.inv - inv0}) "
              f"grid_changed={not np.array_equal(popped._grid, bpd._grid) if popped._grid.shape == bpd._grid.shape else 'shape'}")
    except Exception as e:  # noqa: BLE001
        print(f"    pop_op raised {type(e).__name__}: {e}")
    if db is not None:
        r, c = db
        opts = []
        for t in (TileType.CROSS, TileType.VERT, TileType.HORIZ, TileType.ELBOW_NW, TileType.ELBOW_SE):
            new = bpd.copy()
            new._grid[r, c] = t
            new.rebuild()
            if new.is_valid:
                opts.append((t.name, _num_blanks(new), new.perm.trimcode, new.perm.inv))
        print(f"    valid subs at {db}: {opts}")


def un_gao_huang_absorb(bpd: BPD, cap: int = 60, verbose: bool = False):
    """Blank-driven pop loop with absorption detection.

    Returns ``(letters, status)`` where ``letters`` is a list of
    ``(reflection, row, kind)`` with ``kind in {"grow", "absorb"}`` in the order
    they were removed (i.e. reverse of insertion).
    """
    work = bpd.resize(max(bpd.rows, len(bpd.perm)))
    letters: list[tuple[int, int, str]] = []
    seen: set = set()
    steps = 0

    while _num_blanks(work) > 0 and steps < cap:
        steps += 1
        inv_before = work.perm.inv
        key = _grid_key(work)

        if inv_before > 0 and key not in seen:
            # try a genuine inversion-reducing pop
            seen.add(key)
            popped, (refl, row) = work.pop_op()
            if popped.perm.inv == inv_before - 1:
                letters.append((int(refl), int(row), "grow"))
                work = popped
                seen = set()  # progress made, reset cycle memory
                continue
            # pop failed to reduce inv -> fall through to absorption

        # absorption: the driving blank is an excess (double) crossing
        db = _driving_blank(work)
        if db is None:
            break
        r, c = db
        # the reflection/row this excess blank encodes, read the way pop_op does:
        # pop_op reports (col+1, row+1) for the final mark; here we tag by the
        # blank's own coordinates as a first approximation.
        neutral = _neutralise_blank(work, r, c)
        if neutral is None:
            return letters, f"STUCK@blank({r},{c})"
        letters.append((int(c + 1), int(r + 1), "absorb"))
        work = neutral
        seen = set()

    if _num_blanks(work) > 0:
        return letters, "CAP"
    letters.reverse()
    return letters, "ok"


def letters_to_rows(letters, nrows):
    rows = [[] for _ in range(nrows)]
    for refl, row, _kind in reversed(letters):
        rows[row - 1] = rows[row - 1] + [refl]
    return tuple(tuple(r) for r in rows)


def analyze(codelist):
    w = uncode(codelist)
    n = len(w)
    ell = w.inv
    allwc = WCGraph.all_wc_graphs(w, n)
    wc_rows = {tuple(tuple(r) for r in wc): wc for wc in allwc}
    bpds = BPD.all_unreduced_bpds(w, n)
    unred = [b for b in bpds if not b.is_reduced]
    print(f"\n=== w={w.trimcode} ell={ell}: {len(bpds)} BPDs, {len(unred)} unreduced, {len(allwc)} WCGraphs ===")
    for b in unred:
        blanks = _num_blanks(b)
        diagnose(b)
        letters, status = un_gao_huang_absorb(b)
        ngrow = sum(1 for *_, k in letters if k == "grow")
        nabs = sum(1 for *_, k in letters if k == "absorb")
        rows = letters_to_rows(letters, b.rows) if status == "ok" else None
        is_wc = rows in wc_rows if rows is not None else False
        print(f"  blanks={blanks} status={status:14s} grow={ngrow} absorb={nabs} letters={letters}")
        if rows is not None:
            print(f"      -> rows={rows}  validWC={is_wc}")


if __name__ == "__main__":
    cases = [[1, 0, 2, 1], [2, 0, 2], [1, 0, 3]]
    if len(sys.argv) > 1:
        cases = [[int(c) for c in arg.split(",")] for arg in sys.argv[1:]]
    for c in cases:
        analyze(c)
