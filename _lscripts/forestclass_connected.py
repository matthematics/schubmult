"""Are forest classes of RC graphs connected by adjacent swaps within the class?

Lemma forestlift chains a forest equivalence into single Nadeau-Tewari swaps and
applies the single-swap lemma at each step. That requires every intermediate word
to be the word of an RC graph of the same height, which is not automatic: an
adjacent transposition of a word need not admit a compatible sequence bounded by
that height.

This tests, for each forest class of RC graphs of a fixed height, whether the set
of words occurring in the class is connected under adjacent transpositions that
stay inside the set.
"""

import sys
from collections import Counter, defaultdict

from schubmult import Permutation, RCGraph


def neighbors(word, pool):
    out = []
    for i in range(len(word) - 1):
        if word[i] == word[i + 1]:
            continue
        cand = word[:i] + (word[i + 1], word[i]) + word[i + 2 :]
        if cand in pool:
            out.append(cand)
    return out


def main(n):
    stats = Counter()
    broken = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = max(w.descents()) + 1
        for h in range(md, n + 1):
            classes = defaultdict(set)
            for B in RCGraph.all_rc_graphs(w, h):
                classes[B.forest_invariant].add(B.perm_word)
            for inv, words in classes.items():
                if len(words) <= 1:
                    stats["singleton class"] += 1
                    continue
                start = next(iter(words))
                seen = {start}
                stack = [start]
                while stack:
                    u = stack.pop()
                    for v in neighbors(u, words):
                        if v not in seen:
                            seen.add(v)
                            stack.append(v)
                if seen == words:
                    stats["class connected"] += 1
                else:
                    stats["class DISCONNECTED"] += 1
                    if len(broken) < 5:
                        broken.append((w, h, sorted(words), sorted(seen)))

    print(f"n = {n}")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if broken:
        print("\ndisconnected classes:")
        for w, h, words, reached in broken:
            print(f"  w={w} h={h}")
            print(f"    class words: {words}")
            print(f"    reachable  : {reached}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
