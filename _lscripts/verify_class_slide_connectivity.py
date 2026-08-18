"""Is a forest class connected by SINGLE anti-diagonal (chute/ladder) moves?

If yes, the whole of Lemma forestlift closes as follows:

  (1) forest classes, at fixed permutation u and height m, are connected by
      single anti-diagonal crossing slides that stay inside the class;
  (2) Z intertwines anti-diagonal slides (verified 506/506);
  (3) hence each step of the chain lifts to the same slide upstairs, which
      therefore preserves the forest class upstairs;
  (4) chaining gives B =_forest B', which is Lemma forestlift.

Here we test (1), and also the weaker variant allowing any single-crossing move.
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(perm):
    return len(perm.trimcode)


def rc_all(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h]


def crossings(R):
    out = set()
    for i, row in enumerate(R, start=1):
        for L in row:
            out.add((i, L - i + 1))
    return out


def is_slide(R1, R2):
    """single crossing moved along an anti-diagonal (letter preserved)"""
    c1, c2 = crossings(R1), crossings(R2)
    gone, came = c1 - c2, c2 - c1
    if len(gone) != 1 or len(came) != 1:
        return False
    (i, j), = gone
    (ip, jp), = came
    return i + j == ip + jp


def is_single(R1, R2):
    c1, c2 = crossings(R1), crossings(R2)
    return len(c1 - c2) == 1 and len(c2 - c1) == 1


def connected(nodes, adj):
    if not nodes:
        return True
    start = next(iter(nodes))
    seen = {start}
    stack = [start]
    while stack:
        x = stack.pop()
        for y in adj.get(x, ()):
            if y in nodes and y not in seen:
                seen.add(y)
                stack.append(y)
    return seen == nodes


def main(N=5):
    tot = 0
    bad_slide = 0
    bad_single = 0
    size_hist = Counter()
    examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            u = Permutation(list(pw))
            if maxd(u) == 0:
                continue
            for m in range(max(maxd(u), 1), n + 1):
                if (u, m) in seen:
                    continue
                seen.add((u, m))
                graphs = rc_all(u, m)
                if len(graphs) < 2:
                    continue

                adj_slide = defaultdict(list)
                adj_single = defaultdict(list)
                for a in range(len(graphs)):
                    for b in range(len(graphs)):
                        if a == b:
                            continue
                        R1, R2 = graphs[a], graphs[b]
                        if is_slide(R1, R2):
                            adj_slide[R1].append(R2)
                        if is_single(R1, R2):
                            adj_single[R1].append(R2)

                byclass = defaultdict(set)
                for R in graphs:
                    byclass[R.forest_invariant].add(R)

                for phi, blk in byclass.items():
                    tot += 1
                    size_hist[len(blk)] += 1
                    if not connected(blk, adj_slide):
                        bad_slide += 1
                        if len(examples) < 5:
                            examples.append(
                                (tuple(u[: n]), m, len(blk), [(tuple(rr) for rr in R) and [tuple(r) for r in R] for R in blk])
                            )
                    if not connected(blk, adj_single):
                        bad_single += 1

    print(f"N={N}: forest classes (fixed perm, fixed height) tested = {tot}")
    print(f"  class size histogram = {dict(sorted(size_hist.items()))}")
    print(f"  connected by single ANTI-DIAGONAL slides: {tot - bad_slide}/{tot} (bad {bad_slide})")
    print(f"  connected by any single-crossing move:    {tot - bad_single}/{tot} (bad {bad_single})")
    for e in examples:
        print("   DISCONNECTED:", e)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
