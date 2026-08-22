"""Base-case test with the validated rootlist convention.

NT's insertion order is the reverse of ours, so the prefix U of Proposition 5.8
is the reversal of our word's suffix (the letters after position p+1). The rootlist
is simulated by Lemma 5.7.

Reports: (i) validation of the criterion against forest equivalence; (ii) how
often the bump leaves that suffix untouched, so that the upstairs and downstairs
rootlists coincide; (iii) whether the separating pair survives.
"""

import itertools
import sys
from collections import Counter

from schubmult import Permutation, RCGraph

LO, HI = -4, 40


def injectify(word):
    seen = Counter()
    out = []
    for v in word:
        seen[v] += 1
        out.append((v, seen[v]))
    return out


def rootlist(letters):
    rl = sorted((i, 0) for i in range(LO, HI))
    for a in letters:
        below = [r for r in rl if r < a]
        above = [r for r in rl if r > a]
        rl = sorted(set(rl) - {below[-1], above[0]} | {a})
    return rl


def separators(a, b, rl):
    lo, hi = (a, b) if a < b else (b, a)
    return [r for r in rl if lo < r < hi]


def maxd(perm):
    return max(perm.descents()) + 1 if perm.inv else 0


def main(n):
    stats = Counter()
    broken = []

    for w in Permutation.all_permutations(n):
        if w.inv == 0:
            continue
        md = maxd(w)
        for h in range(md, n + 1):
            if md <= h - 1:
                continue
            graphs = [B for B in RCGraph.all_rc_graphs(w, h) if len(B[-1]) == 0]
            images = {}
            for B in graphs:
                try:
                    images[B] = B.zero_out_last_row()
                except Exception:
                    pass
            for B1, B2 in itertools.permutations(images, 2):
                R1, R2 = images[B1], images[B2]
                if R1.perm != R2.perm:
                    continue
                d1, d2 = R1.perm_word, R2.perm_word
                diff = [i for i in range(len(d1)) if d1[i] != d2[i]]
                if len(diff) != 2 or diff[1] != diff[0] + 1:
                    continue
                p = diff[0]
                if not (d1[p] == d2[p + 1] and d1[p + 1] == d2[p]):
                    continue

                inj = injectify(d1)
                rl_down = rootlist(list(reversed(inj[p + 2 :])))
                seps = separators(inj[p], inj[p + 1], rl_down)
                actual = R1.forest_invariant == R2.forest_invariant
                stats["validation agree" if (len(seps) >= 2) == actual else "validation MISMATCH"] += 1
                if not actual:
                    continue

                up = B1.perm_word
                tail_untouched = all(up[i] == d1[i] for i in range(p + 2, len(up)))
                stats["feq: suffix untouched (base case)" if tail_untouched
                      else "feq: suffix touched"] += 1

                inj_up = injectify(up)
                rl_up = rootlist(list(reversed(inj_up[p + 2 :])))
                seps_up = separators(inj_up[p], inj_up[p + 1], rl_up)
                tag = "base" if tail_untouched else "step"
                stats[f"  {tag}: upstairs separated" if len(seps_up) >= 2
                      else f"  {tag}: upstairs NOT separated"] += 1
                if tail_untouched:
                    stats[f"  base: rootlists equal = {rl_up == rl_down}"] += 1
                    stats[f"  base: original separators surviving = {len([r for r in seps if r in seps_up])}"] += 1
                if len(seps_up) < 2 and len(broken) < 5:
                    broken.append((up, d1, p, seps, seps_up))

    print(f"n = {n}")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if broken:
        print("\nupstairs-not-separated cases:")
        for up, d1, p, s, su in broken:
            print(f"  up={up} down={d1} p={p} down_seps={s} up_seps={su}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
