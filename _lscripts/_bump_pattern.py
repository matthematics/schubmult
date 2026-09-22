"""Pattern of a downward Little bump at position with right root (a,b): express w^{-1} w' as transpositions."""
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

def perm_arr(word, n):
    w = Permutation.ref_product(*word) if word else Permutation([])
    arr = list(w); arr += list(range(len(arr) + 1, n + 1))
    return arr[:n]

def reduced_words(w):
    if w.inv == 0:
        yield ()
        return
    for d in w.descents():
        sw = w * Permutation.ref_product(d + 1)
        for rw in reduced_words(sw):
            yield (*rw, d + 1)

random.seed(0)
pat = Counter()
n = 9
for N in range(3, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 7:
            continue
        words = list(reduced_words(w))
        if len(words) > 30:
            words = random.sample(words, 30)
        for word in words:
            for q in range(len(word)):
                a, b = rroot(word, q)
                nw = ltb_down(word, q)
                if nw is None:
                    continue
                wa, wb = perm_arr(word, n), perm_arr(nw, n)
                diff = [i + 1 for i in range(n) if wa[i] != wb[i]]
                # classify
                if len(diff) == 2:
                    pat[("2 positions", "is (a,b)" if tuple(diff) == (a, b) else "other")] += 1
                elif len(diff) == 3:
                    c = [x for x in diff if x not in (a, b)]
                    if len(c) == 1 and a in diff and b in diff:
                        c = c[0]
                        # which of a,b keeps its value?
                        keep = "a fixed" if wa[a-1] == wb[a-1] else "b fixed" if wa[b-1] == wb[b-1] else "neither"
                        pat[("3 positions", "c<a" if c < a else "a<c<b" if c < b else "c>b", keep)] += 1
                    else:
                        pat[("3 positions", "not containing a,b")] += 1
                else:
                    pat[(len(diff), "positions")] += 1
print(pat)
