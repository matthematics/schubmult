"""Check the run bookkeeping in the proof of T1: Little bump of word(R) at the last position (x = leftmost cross of
last nonempty row alpha, column j >= 2).  Each maximal run inside the prefix starts at a position k whose root in the
current prefix is {a-1, a} where a is the last letter before the decrement; check min root >= alpha, and that the
prefix restricted runs are complete Little bumps of the prefix (prefix reduced at run end)."""
import sys
from collections import Counter
from itertools import permutations

from _lb_common import rroots
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import is_reduced

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
stats = Counter()
for N0 in range(2, NMAX + 1):
    for a_ in permutations(range(1, N0 + 1)):
        w = Permutation(list(a_))
        if w.inv == 0:
            continue
        for R in RCGraph.all_rc_graphs(w, N0):
            rows = [tuple(r) for r in R]
            nonempty = [i for i, r in enumerate(rows) if r]
            alpha = nonempty[-1] + 1
            row = rows[alpha - 1]
            j = row[-1] - alpha + 1
            if j < 2:
                continue
            word = [l for r in rows for l in r]
            last = len(word) - 1
            cur = last
            in_prefix = False
            while True:
                if cur == last:
                    a = word[cur]
                    word[cur] -= 1
                    if is_reduced(word):
                        break
                    rts = rroots(word)
                    partners = [q for q in range(len(word)) if q != cur and set(rts[q]) == set(rts[cur])]
                    assert len(partners) == 1
                    k = partners[0]
                    # root of k in the prefix alone
                    pref = word[:last]
                    prts = rroots(pref)
                    s, beta = sorted(prts[k])
                    stats[("run root == {a-1,a}", {s, beta} == {a - 1, a})] += 1
                    stats[("s >= alpha", s >= alpha)] += 1
                    stats[("prefix deletable at k", is_reduced(pref[:k] + pref[k + 1 :]))] += 1
                    cur = k
                else:
                    word[cur] -= 1
                    assert word[cur] >= 1
                    if is_reduced(word):
                        stats["prefix reduced at end of run"] += 1 if is_reduced(word[:last]) else 0
                        break
                    rts = rroots(word)
                    partners = [q for q in range(len(word)) if q != cur and set(rts[q]) == set(rts[cur])]
                    assert len(partners) == 1
                    nxt = partners[0]
                    if nxt == last:
                        stats[("prefix reduced when run leaves", is_reduced(word[:last]))] += 1
                    cur = nxt
            stats["cases"] += 1
for k, v in sorted(stats.items(), key=str):
    print(v, k)
