"""Interaction of two p-high bumps: chain overlap vs commutation."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def high_bumps(word, p):
    out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            res = ltb_down(word, q, return_chain=True)
            if res[0] is not None:
                out[(a, b)] = (res[0], res[1], q)
    return out

random.seed(0)
stats = Counter()
ex = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 15:
            words = random.sample(words, 15)
        for word in words:
            for p in range(1, maxd(w)):
                hb = high_bumps(word, p)
                keys = sorted(hb)
                for i in range(len(keys)):
                    for j in range(i + 1, len(keys)):
                        b1, b2 = keys[i], keys[j]
                        w1, ch1, q1 = hb[b1]; w2, ch2, q2 = hb[b2]
                        overlap = bool(set(ch1) & set(ch2))
                        h1 = high_bumps(w1, p); h2 = high_bumps(w2, p)
                        # does b2 still exist in w1 as a p-high inversion, and does bumping give same as b1 on w2?
                        commute = (b2 in h1 and b1 in h2 and h1[b2][0] == h2[b1][0])
                        # alternative join: some single bumps join
                        join1 = bool({x[0] for x in h1.values()} & {x[0] for x in h2.values()}) or w1 == w2
                        stats[("overlap" if overlap else "disjoint", "commute" if commute else "nocommute", "join1" if join1 else "nojoin1")] += 1
                        if not overlap and not commute and len(ex) < 6:
                            ex.append((word, p, b1, ch1, b2, ch2, w1, w2))
print(stats)
for e in ex:
    print(e)
