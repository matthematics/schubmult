"""Confluence of non-degenerate p-high downward Little bumps (word level), and stats on degenerate cases."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def high_bumps(word, p):
    out = set(); degen = 0
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is None:
                degen += 1
            else:
                out.add(nw)
    return out, degen

memo = {}
def normal_forms(word, p):
    key = (word, p)
    if key in memo:
        return memo[key]
    nb, _ = high_bumps(word, p)
    res = frozenset([word]) if not nb else frozenset().union(*(normal_forms(x, p) for x in nb))
    memo[key] = res
    return res

random.seed(0)
stats = Counter(); degs = Counter(); nfmax = Counter()
bad = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 25:
            words = random.sample(words, 25)
        for word in words:
            for p in range(1, maxd(w)):
                _, d = high_bumps(word, p)
                degs[d > 0] += 1
                nf = normal_forms(word, p)
                stats[len(nf)] += 1
                for x in nf:
                    nfmax[maxd(Permutation.ref_product(*x)) <= p] += 1
                if len(nf) > 1 and len(bad) < 5:
                    bad.append((word, p, sorted(nf)))
print("normal forms:", stats)
print("words having degenerate high bumps:", degs)
print("normal forms with maxd<=p:", nfmax)
for b in bad:
    print(b)
