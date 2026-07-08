"""Calibrate a BPD-resize-based WCGraph.zero_out_last_row against the existing
implementation, using the MBPD<->WCGraph bijection.

Pipeline:  WCGraph --Psi--> MBPD --tiles--> BPD --resize(len(wc)-1)--> BPD
           --tiles--> MBPD --Phi--> WCGraph  (then normalize trailing empties).
"""

import itertools

import numpy as np

from schubmult.combinatorics.bpd import BPD, TileType
from schubmult.combinatorics.mbpd import MBPD
from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.wc_graph import WCGraph

M2B = {
    "B": TileType.BLANK, "P": TileType.CROSS, "H": TileType.HORIZ,
    "V": TileType.VERT, "J": TileType.ELBOW_NW, "R": TileType.ELBOW_SE,
    "M": TileType.ELBOW_NW,
}
B2M = {
    TileType.BLANK: "B", TileType.CROSS: "P", TileType.HORIZ: "H",
    TileType.VERT: "V", TileType.ELBOW_NW: "J", TileType.ELBOW_SE: "R",
}


def mbpd_to_bpd(D):
    g = np.array([[M2B[D.tile(i, j)] for j in range(1, D.n + 1)] for i in range(1, D.n + 1)], dtype=TileType)
    marks = frozenset((i, j) for i in range(1, D.n + 1) for j in range(1, D.n + 1) if D.tile(i, j) == "M")
    return BPD(g), marks


def bpd_to_mbpd(b, marks):
    N = max(b.rows, len(b.perm))
    b = b.resize(N)
    grid = [[B2M[b._grid[i, j]] for j in range(N)] for i in range(N)]
    D = MBPD.from_tiles(grid)
    for (i, j) in marks:
        if i <= N and j <= N and D.tile(i, j) == "J":
            D = D.with_tile(i, j, "M")
    return D


def new_zero(wc):
    R = len(wc) - 1
    D = wc.to_mbpd()
    b, marks = mbpd_to_bpd(D)
    bz = b.resize(R)
    Dz = bpd_to_mbpd(bz, frozenset((i, j) for (i, j) in marks if i <= R))
    return Dz.phi().to_wcgraph().resize(R)


def norm(wc):
    t = [tuple(r) for r in wc]
    while t and t[-1] == ():
        t.pop()
    return tuple(t)


def main():
    tot = same = diff = total_new = none_old = 0
    wt_fail = perm_ok = 0
    for n in (3, 4, 5):
        for pl in itertools.permutations(range(1, n + 1)):
            w = Permutation(list(pl))
            for wc in WCGraph.all_wc_graphs(w):
                if len(wc) == 0 or len(wc[-1]) != 0:
                    continue
                total_new += 1
                new = new_zero(wc)
                # weight of a WCGraph = multiset of compatible values; zero_out
                # drops the (empty) last row so total weight is preserved.
                if not new.is_valid:
                    print(f"INVALID new for wc={norm(wc)} -> {norm(new)}")
                if sum(new.weight) == sum(wc.weight):
                    perm_ok += 1
                else:
                    wt_fail += 1
                    if wt_fail <= 10:
                        print(f"WEIGHT wc={norm(wc)} wt={wc.weight} -> new={norm(new)} wt={new.weight}")
                old = wc.zero_out_last_row()
                if old is None:
                    none_old += 1
                    continue
                tot += 1
                if norm(new) == norm(old):
                    same += 1
                else:
                    diff += 1
                    if diff <= 15:
                        print(f"DIFF pl={pl} wc={norm(wc)} old={norm(old)} oldperm={list(old.perm)} new={norm(new)} newperm={list(new.perm)}")
    print(f"\ntotal empty-last-row WCGraphs={total_new}")
    print(f"weight preserved: {perm_ok}  (weight changed: {wt_fail})")
    print(f"old!=None: {tot}  (same={same}, diff={diff})")
    print(f"old==None (non-core-reduced, new is TOTAL): {none_old}")


if __name__ == "__main__":
    main()
