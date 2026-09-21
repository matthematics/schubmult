"""Is the relative order pattern of the word letters (sign(a_i - a_j) for all i<j) preserved by p-high bumps?
Also test for Little bumps in general (any root)."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def sgn(x):
    return (x > 0) - (x < 0)

def pattern(word):
    n = len(word)
    return tuple(sgn(word[i] - word[j]) for i in range(n) for j in range(i + 1, n))

random.seed(0)
stats = Counter(); ex = []
for N in range(3, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 8:
            continue
        words = list(reduced_words(w))
        if len(words) > 15:
            words = random.sample(words, 15)
        for word in words:
            pat = pattern(word)
            for q in range(len(word)):
                a, b = rroot(word, q)
                nw = ltb_down(word, q)
                if nw is None:
                    continue
                same = pattern(nw) == pat
                for p in range(1, maxd(w)):
                    stats[("p-high" if a > p else "not p-high", same)] += 1
                if not same and a > 1 and len(ex) < 5:
                    ex.append((word, q, (a, b), nw))
print(stats)
for e in ex:
    print(e)
