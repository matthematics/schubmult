"""Separate the WORD level from the RC GRAPH level.

Forest equivalence is defined on WORDS.  Nadeau--Tewari say the forest class of a
word is connected under the permitted local swaps.  A forest class of RC GRAPHS
is bigger: several graphs can share a word (different compatible sequences), and
no word swap connects those.

So the reduction to the simple case has two halves:

  (W) words: connect word(R1) to word(R2) by permitted local swaps inside the
      class -- this is Nadeau--Tewari;
  (G) graphs: for graphs sharing a word, something else is needed.

But (G) is FREE for forestlift: if word(R1) == word(R2) then R1 =_forest R2
holds trivially, and what we need upstairs is word(B1) =_forest word(B2).

This script measures the true scope of forestlift on the nontrivial branch:
  (1) all pairs B, B' in RC_h(w)^0 with Z(B) =_forest Z(B'): does B =_forest B'?
  (2) of those, how many have word(Z(B)) and word(Z(B')) differing by a single
      swap (the simple case), and how many need a longer chain?
  (3) within the class, is the set of words of Z-images connected by single
      swaps?
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


def swaps(w1, w2):
    if len(w1) != len(w2):
        return []
    return [p for p in range(len(w1) - 1)
            if w1[:p] + (w1[p + 1], w1[p]) + w1[p + 2:] == w2 and abs(w1[p] - w1[p + 1]) >= 2]


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
    lift_tot = lift_ok = 0
    dist = Counter()
    simple_tot = simple_ok = 0
    sameword_tot = sameword_ok = 0
    wconn_tot = wconn_ok = 0
    wclass_size = Counter()
    lastpos = Counter()
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

            # (1) forestlift, full scope
            for R1 in Rs:
                for R2 in Rs:
                    if R1 is R2:
                        continue
                    if R1.forest_invariant != R2.forest_invariant:
                        continue
                    B1, B2 = Zinv[R1], Zinv[R2]
                    lift_tot += 1
                    ok = B1.forest_invariant == B2.forest_invariant
                    if ok:
                        lift_ok += 1
                    w1, w2 = word(R1), word(R2)
                    sp = swaps(w1, w2)
                    if w1 == w2:
                        dist["same word"] += 1
                        sameword_tot += 1
                        if ok:
                            sameword_ok += 1
                    elif sp:
                        dist["one swap"] += 1
                        simple_tot += 1
                        if ok:
                            simple_ok += 1
                        for p in sp:
                            lastpos[p == len(w1) - 2] += 1
                    else:
                        dist["longer"] += 1
                    if not ok and len(examples) < 5:
                        examples.append((tuple(w[:n]), h, w1, w2,
                                         [tuple(r) for r in B1], [tuple(r) for r in B2]))

            # (3) word-level connectivity inside a class
            byclass = defaultdict(set)
            for R in Rs:
                byclass[(R.perm, R.forest_invariant)].add(word(R))
            for key, ws in byclass.items():
                wconn_tot += 1
                wclass_size[len(ws)] += 1
                adj = defaultdict(list)
                for a in ws:
                    for b in ws:
                        if a != b and swaps(a, b):
                            adj[a].append(b)
                if connected(ws, adj):
                    wconn_ok += 1

    print(f"N={N}  (h = maxd(w), nontrivial branch)")
    print(f"  FORESTLIFT, all pairs: {lift_ok}/{lift_tot}")
    print(f"     breakdown of the downstairs relation: {dict(dist)}")
    print(f"       same word:  {sameword_ok}/{sameword_tot}")
    print(f"       one swap:   {simple_ok}/{simple_tot}   (swap at last two letters? {dict(lastpos)})")
    print()
    print(f"  WORD classes (of Z-images) connected by single swaps: {wconn_ok}/{wconn_tot}")
    print(f"     distinct words per class: {dict(sorted(wclass_size.items()))}")
    for e in examples:
        print("\n   FORESTLIFT FAILS:")
        print(f"     w={e[0]} h={e[1]}  words {e[2]} / {e[3]}")
        print(f"     B1={e[4]}  B2={e[5]}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
