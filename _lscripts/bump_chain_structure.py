"""Structure of the Little bump chain, at the level of words.

The bump decrements the letter at some position; if the result is not reduced
there is a unique other position whose deletion restores reducedness, and the
bump continues there. This records the chain of positions and tests two
predictions of the root computation:

  collision at j requires   M^{-1}(alpha_{a_j}) = -alpha_{l-1},

where l is the value being decremented and M is the product of the letters lying
strictly between j and the bumped position. Taking j = p-1 makes M empty and the
condition reads alpha_{a_j} = -alpha_{l-1}, which is impossible for positive
roots; so the chain should never step to the immediately preceding position.
"""

import sys
from collections import Counter

from schubmult import Permutation, RCGraph


def perm_of(word, n):
    p = list(range(n + 2))
    for a in word:
        p[a - 1], p[a] = p[a], p[a - 1]
    return tuple(p)


def inv_count(p):
    return sum(1 for i in range(len(p)) for j in range(i + 1, len(p)) if p[i] > p[j])


def is_reduced(word, n):
    return len(word) == inv_count(perm_of(word, n))


def bump_chain(word, start, n):
    """Little bump beginning at `start`; returns the list of positions decremented."""
    w = list(word)
    chain = [start]
    i = start
    for _ in range(200):
        if w[i] == 1:
            w = [x + 1 if k != i else x for k, x in enumerate(w)]
        else:
            w[i] -= 1
        if is_reduced(w, n):
            return chain, tuple(w)
        cands = [j for j in range(len(w)) if j != i and is_reduced(w[:j] + w[j + 1 :], n)]
        if len(cands) != 1:
            return None, None
        i = cands[0]
        chain.append(i)
    return None, None


def main(n):
    stats = Counter()
    rightward = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = max(w.descents()) + 1
        for h in range(md, n + 1):
            if md != h:
                continue
            for B in RCGraph.all_rc_graphs(w, h):
                if len(B[-1]) != 0:
                    continue
                word = B.perm_word
                # the bump begins at the letter whose root is the last descent
                roots = [B.left_to_right_inversion(i) for i in range(B.perm.inv)]
                try:
                    start = roots.index((h, h + 1))
                except ValueError:
                    continue
                chain, out = bump_chain(word, start, n + 3)
                if chain is None:
                    stats["chain failed"] += 1
                    continue
                stats[f"chain length {len(chain)}"] += 1
                for u, v in zip(chain, chain[1:]):
                    stats["step leftward" if v < u else "step RIGHTWARD"] += 1
                    if v == u - 1:
                        stats["step to IMMEDIATE predecessor"] += 1
                    if v > u and len(rightward) < 5:
                        rightward.append((word, chain))

    print(f"n = {n}")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if rightward:
        print("\nrightward steps:")
        for word, chain in rightward:
            print(f"  {word}: chain {chain}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
