"""Window lemma restricted to the configurations that actually arise.

Only the suffix U (reversed, in NT insertion order) matters for the separation
criterion, and only the decrements the bump performs inside that suffix. Tests
whether rl(U') and rl(U) agree outside the windows [v-1, v+1] of the values
decremented within the suffix, over the forest-equivalent instances.
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
    escapes = []

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
                if R1.forest_invariant != R2.forest_invariant:
                    continue

                up, down = B1.perm_word, d1
                inj_dn, inj_up = injectify(down), injectify(up)
                U_dn = list(reversed(inj_dn[p + 2 :]))
                U_up = list(reversed(inj_up[p + 2 :]))
                rl_dn, rl_up = rootlist(U_dn), rootlist(U_up)

                dec_vals = {down[i] for i in range(p + 2, len(down)) if up[i] != down[i]}
                windows = set()
                for v in dec_vals:
                    windows |= {v - 1, v, v + 1}

                sym = set(rl_dn) ^ set(rl_up)
                outside = sorted(r for r in sym if r[0] not in windows)
                stats[f"suffix decrements = {len(dec_vals)}"] += 1
                if outside:
                    stats["ESCAPES windows"] += 1
                    if len(escapes) < 5:
                        escapes.append((up, down, p, sorted(dec_vals), outside))
                else:
                    stats["agrees outside windows"] += 1

                # the count invariant
                s_dn = len(separators(inj_dn[p], inj_dn[p + 1], rl_dn))
                s_up = len(separators(inj_up[p], inj_up[p + 1], rl_up))
                stats[f"separator count {s_dn} -> {s_up}"] += 1
                if s_up < s_dn:
                    stats["count DECREASED"] += 1

    print(f"n = {n}  (forest-equivalent instances)")
    for k in sorted(stats):
        print(f"  {k}: {stats[k]}")
    if escapes:
        print("\nescapes:")
        for up, down, p, vals, outside in escapes:
            print(f"  up={up} down={down} p={p} decremented in suffix {vals}")
            print(f"    outside: {outside}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
