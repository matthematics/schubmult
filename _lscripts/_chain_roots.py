"""At each step of a p-high bump chain, record the right-root (as positions) of the letter about to be
decremented. Test whether both positions are always > p (i.e. p-high bumps never touch crossings involving a wire <= p)."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroots, reduced_words, maxd
from schubmult.utils.perm_utils import is_reduced
from collections import Counter

def chain_roots(word, index):
    """Downward chain; returns list of (root_before_decrement) or None if degenerate."""
    word = [*word]
    out = []
    while True:
        if word[index] == 1:
            return None
        rts = rroots(word)
        a, b = rts[index]
        out.append((min(a, b), max(a, b)))
        word[index] -= 1
        if is_reduced(word):
            return out, tuple(word)
        rts = rroots(word)
        target = set(rts[index])
        partners = [i for i in range(len(word)) if i != index and set(rts[i]) == target]
        if len(partners) != 1:
            return None
        index = partners[0]

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
            rts = rroots(word)
            for p in range(1, maxd(w)):
                for q in range(len(word)):
                    a, b = sorted(rts[q])
                    if a <= p:
                        continue
                    res = chain_roots(word, q)
                    if res is None:
                        continue
                    roots, nw = res
                    allhigh = all(r[0] > p for r in roots)
                    # also: does the chain keep the smaller wire a fixed?
                    keep_a = all(a in r for r in roots)
                    keep_some = all(len(set(roots[0]) & set(r)) > 0 for r in roots)
                    stats[("all roots p-high", allhigh)] += 1
                    stats[("a in every root", keep_a)] += 1
                    stats[("shares wire with first", keep_some)] += 1
                    if not allhigh and len(ex) < 5:
                        ex.append((word, p, q, roots))
print(stats)
for e in ex:
    print(e)
