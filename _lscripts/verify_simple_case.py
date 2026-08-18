"""The simple case: a SINGLE local swap, anywhere in the word.

The user's point: a forest class is a connected graph under the generating local
swaps, so Lemma forestlift only needs the simple case -- one swap -- and the rest
follows by chaining along the graph.  No reduction to the last two letters, and
no interior antidiagonal move apparatus, is required.

SIMPLE CASE.  Let R1, R2 be bounded RC graphs of the same height m and the same
permutation, with word(R2) obtained from word(R1) by swapping two adjacent
commuting letters, and R1 =_forest R2.  Let B1, B2 be their preimages under Z in
RC_{m+1}(w)^0 for a common w.  Then B1 =_forest B2.

Chaining this along a connecting path inside the class gives forestlift.

Tested on the nontrivial branch h = maxd(w), where Algorithm 3 really runs.
We also record:
  - connectivity of each class under single swaps that stay in the class,
    restricted to the graphs that actually occur as Z-images;
  - whether the swap positions needed are ever forced to be the last two.
"""

import sys
from collections import Counter, defaultdict
from itertools import permutations

from schubmult import Permutation
from schubmult.combinatorics.rc_graph import RCGraph


def maxd(p):
    return len(p.trimcode)


def rc0(perm, h):
    return [R for R in RCGraph.all_rc_graphs(perm, h) if len(R) == h and len(R[-1]) == 0]


def word(R):
    return tuple(R.perm_word)


def swap_positions(w1, w2):
    """positions p where w2 = w1 with entries p,p+1 transposed; else []"""
    if len(w1) != len(w2):
        return []
    out = []
    for p in range(len(w1) - 1):
        cand = w1[:p] + (w1[p + 1], w1[p]) + w1[p + 2 :]
        if cand == w2 and abs(w1[p] - w1[p + 1]) >= 2:
            out.append(p)
    return out


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


def main(N=6):
    simple_tot = simple_ok = 0
    pos_hist = Counter()
    lasttwo = Counter()
    conn_tot = conn_ok = 0
    classsize = Counter()
    examples = []

    seen = set()
    for n in range(2, N + 1):
        for pw in permutations(range(1, n + 1)):
            w = Permutation(list(pw))
            h = maxd(w)
            if h < 2 or (w, h) in seen:
                continue
            seen.add((w, h))

            graphs = rc0(w, h)
            if len(graphs) < 2:
                continue
            Zinv = {}
            for B in graphs:
                Zinv[B.zero_out_last_row()] = B
            Rs = list(Zinv)

            # adjacency within a forest class by single swaps
            adj = defaultdict(list)
            for R1 in Rs:
                for R2 in Rs:
                    if R1 is R2 or R1.perm != R2.perm:
                        continue
                    if R1.forest_invariant != R2.forest_invariant:
                        continue
                    ps = swap_positions(word(R1), word(R2))
                    if not ps:
                        continue
                    adj[R1].append(R2)

                    # THE SIMPLE CASE
                    simple_tot += 1
                    B1, B2 = Zinv[R1], Zinv[R2]
                    if B1.forest_invariant == B2.forest_invariant:
                        simple_ok += 1
                    elif len(examples) < 5:
                        examples.append(
                            (tuple(w[:n]), h, word(R1), word(R2), ps,
                             [tuple(r) for r in B1], [tuple(r) for r in B2])
                        )
                    for p in ps:
                        pos_hist[p - (len(word(R1)) - 2)] += 1
                        lasttwo[p == len(word(R1)) - 2] += 1

            # connectivity of each class, among the graphs occurring as Z-images
            byclass = defaultdict(set)
            for R in Rs:
                byclass[(R.perm, R.forest_invariant)].add(R)
            for key, blk in byclass.items():
                conn_tot += 1
                classsize[len(blk)] += 1
                if connected(blk, adj):
                    conn_ok += 1

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  SIMPLE CASE (one swap, class preserved) => upstairs forest equivalent: {simple_ok}/{simple_tot}")
    print()
    print(f"  classes connected by single in-class swaps: {conn_ok}/{conn_tot}")
    print(f"     class size histogram: {dict(sorted(classsize.items()))}")
    print()
    print(f"  swap position relative to the END of the word (0 = last two): {dict(sorted(pos_hist.items()))}")
    print(f"     is the swap at the last two letters? {dict(lasttwo)}")
    for e in examples:
        print("\n   SIMPLE CASE FAILS:")
        print(f"     w={e[0]} h={e[1]}  {e[2]} -> {e[3]}  at positions {e[4]}")
        print(f"     B1={e[5]}  B2={e[6]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
