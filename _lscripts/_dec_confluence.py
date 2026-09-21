"""Confluence of the sub-system of p-high bumps with strictly decreasing chains (no stack pushes), word level.
Also: within this sub-system, is local confluence always a *commuting diamond* (join depth 1)?"""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def dec_bumps(word, p):
    out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw, ch = ltb_down(word, q, return_chain=True)
            if nw is None:
                continue
            if all(ch[k] > ch[k + 1] for k in range(len(ch) - 1)):
                out[(a, b)] = nw
    return out

memo = {}
def NF(word, p):
    key = (word, p)
    if key in memo:
        return memo[key]
    nb = dec_bumps(word, p)
    res = frozenset([word]) if not nb else frozenset().union(*(NF(x, p) for x in nb.values()))
    memo[key] = res
    return res

random.seed(0)
stats = Counter(); joins = Counter()
bad = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 20:
            words = random.sample(words, 20)
        for word in words:
            for p in range(1, maxd(w)):
                nf = NF(word, p)
                stats[len(nf)] += 1
                if len(nf) > 1 and len(bad) < 5:
                    bad.append((word, p, sorted(nf)))
                B = dec_bumps(word, p)
                keys = sorted(B)
                for i in range(len(keys)):
                    for j in range(i + 1, len(keys)):
                        w1, w2 = B[keys[i]], B[keys[j]]
                        B1, B2 = dec_bumps(w1, p), dec_bumps(w2, p)
                        R1 = set(B1.values()) | {w1}; R2 = set(B2.values()) | {w2}
                        if R1 & R2:
                            comm = keys[j] in B1 and keys[i] in B2 and B1[keys[j]] == B2[keys[i]]
                            joins["commute" if comm else "depth1-other"] += 1
                        else:
                            joins["depth>=2"] += 1
print("NF counts:", stats)
print("joins:", joins)
for b in bad:
    print(b)
