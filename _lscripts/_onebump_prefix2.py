"""Per-bump prefix analysis. For Y (empty last row, maxd = height m) and p < m, let u = word of Y_{<=p},
u' = word of (lmap Y)_{<=p}. Check: u' == u, or u' == ltb_q(u) (word-level downward Little bump) where q is the
largest changed position. Also record the root (a,b) of position q in u and whether a > p."""
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter

def maxd(w):
    return len(w.trimcode)

def ltb_down(word, index):
    word = [*word]
    while True:
        assert word[index] > 1
        word[index] -= 1
        if is_reduced(word):
            return tuple(word)
        index = find_reduced_fail(word, index)

def top_word(X, p):
    rows = [tuple(r) for r in X][:p]
    return tuple(x for r in rows for x in r)

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w)
        if m < 3:
            continue
        for Y in RCGraph.all_rc_graphs(w, m):
            rows = [tuple(r) for r in Y]
            if len(rows[-1]) != 0:
                continue
            Y1 = Y.little_bump_desc()
            for p in range(2, m):
                u = top_word(Y, p)
                u1 = top_word(Y1, p)
                if u == u1:
                    stats["unchanged"] += 1
                    continue
                changed = [k for k in range(len(u)) if u[k] != u1[k]]
                if any(u1[k] != u[k] - 1 for k in changed):
                    stats["FAIL: not -1 changes"] += 1
                    continue
                q = max(changed)
                ok = ltb_down(u, q) == u1
                # root of position q in u (left-to-right)
                T = RCGraph([tuple(r) for r in Y][:p])
                a, b = T.left_to_right_inversion(q)
                stats[("bump-at-entry" if ok else "FAIL", "a>p" if a > p else "a<=p")] += 1
                if not ok and len(bad) < 5:
                    bad.append((rows, p, u, u1, q, (a, b)))
print(stats)
for b in bad:
    print(b)
