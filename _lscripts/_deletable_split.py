"""If x, y reduced words, xy = (reduced word with one letter decremented) non-reduced, can both deletable
positions lie in x?"""
import random
from itertools import permutations
from schubmult import Permutation
from schubmult.utils.perm_utils import is_reduced
from collections import Counter

def reduced_words(w):
    if w.inv == 0:
        yield ()
        return
    for d in w.descents():
        sw = w * Permutation.ref_product(d + 1)
        for rw in reduced_words(sw):
            yield (*rw, d + 1)

random.seed(0)
stats = Counter()
ex = []
for N in range(3, 7):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv == 0 or w.inv > 7:
            continue
        words = list(reduced_words(w))
        if len(words) > 30:
            words = random.sample(words, 30)
        for word in words:
            for k in range(len(word)):
                if word[k] == 1:
                    continue
                nw = [*word]; nw[k] -= 1
                if is_reduced(nw):
                    continue
                dele = [j for j in range(len(nw)) if is_reduced(nw[:j] + nw[j + 1:])]
                assert len(dele) == 2 and k in dele, (word, k, dele)
                for cut in range(1, len(nw)):
                    x, y = nw[:cut], nw[cut:]
                    if not (is_reduced(x) and is_reduced(y)):
                        continue
                    where = tuple(sorted("x" if d < cut else "y" for d in dele))
                    stats[where] += 1
                    if where == ("x", "x") and len(ex) < 5:
                        ex.append((word, k, cut, dele))
print(stats)
for e in ex:
    print(e)
