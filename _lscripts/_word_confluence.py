"""Word-level confluence test: p-high downward Little bumps (at letters whose right-transported root (a,b) has a>p),
disallowing bumps that would create letter 0. Are normal forms unique?"""
import random
from itertools import permutations
from schubmult import Permutation
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter

def rroot(word, q):
    a, b = word[q], word[q] + 1
    for x in word[q + 1:]:
        a = a + 1 if a == x else a - 1 if a == x + 1 else a
        b = b + 1 if b == x else b - 1 if b == x + 1 else b
    return (min(a, b), max(a, b))

def ltb_down(word, index):
    word = [*word]
    while True:
        if word[index] == 1:
            return None
        word[index] -= 1
        if is_reduced(word):
            return tuple(word)
        index = find_reduced_fail(word, index)
        if index is None:
            return None

def high_bumps(word, p):
    out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None:
                out.add(nw)
    return out

def normal_forms(word, p, memo):
    key = (word, p)
    if key in memo:
        return memo[key]
    nb = high_bumps(word, p)
    if not nb:
        res = frozenset([word])
    else:
        res = frozenset().union(*(normal_forms(x, p, memo) for x in nb))
    memo[key] = res
    return res

def reduced_words(w):
    if w.inv == 0:
        yield ()
        return
    for d in w.descents():
        sw = w * Permutation.ref_product(d + 1)
        for rw in reduced_words(sw):
            yield (*rw, d + 1)

def maxd(w):
    return len(w.trimcode)

random.seed(0)
stats = Counter()
bad = []
memo = {}
for N in range(3, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 7:
            continue
        words = list(reduced_words(w))
        if len(words) > 30:
            words = random.sample(words, 30)
        for word in words:
            for p in range(1, maxd(w)):
                nf = normal_forms(word, p, memo)
                stats[len(nf)] += 1
                if len(nf) > 1 and len(bad) < 5:
                    bad.append((word, p, sorted(nf)))
print(stats)
for b in bad:
    print(b)
