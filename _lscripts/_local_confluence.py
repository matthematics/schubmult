"""Local confluence structure of p-high bumps. For each word and each pair of distinct p-high bumps b1,b2,
find minimal joining: depth-1 (diamond), depth-2, etc., within p-high bumps. Also test plain commutation
when position-supports are disjoint."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroot, ltb_down, reduced_words, maxd
from collections import Counter

def perm_arr(word, n=10):
    w = Permutation.ref_product(*word) if word else Permutation([])
    arr = list(w); arr += list(range(len(arr) + 1, n + 1))
    return tuple(arr[:n])

def high_bumps(word, p):
    out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None:
                out[(a, b)] = nw
    return out

def reach(word, p, depth):
    """Set of words reachable by <= depth p-high bumps."""
    cur = {word}; allr = {word}
    for _ in range(depth):
        nxt = set()
        for x in cur:
            nxt |= set(high_bumps(x, p).values())
        allr |= nxt; cur = nxt
    return allr

random.seed(0)
stats = Counter(); comm = Counter()
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
            wa = perm_arr(word)
            for p in range(1, maxd(w)):
                hb = high_bumps(word, p)
                keys = sorted(hb)
                for i in range(len(keys)):
                    for j in range(i + 1, len(keys)):
                        b1, b2 = keys[i], keys[j]
                        w1, w2 = hb[b1], hb[b2]
                        # supports
                        s1 = {k + 1 for k in range(10) if wa[k] != perm_arr(w1)[k]}
                        s2 = {k + 1 for k in range(10) if wa[k] != perm_arr(w2)[k]}
                        disjoint = not (s1 & s2)
                        # joinability depth
                        found = None
                        for d in range(0, 4):
                            if reach(w1, p, d) & reach(w2, p, d):
                                found = d; break
                        stats[("join depth", found)] += 1
                        if disjoint:
                            # commutation: b2 applied to w1 equals b1 applied to w2 ?
                            h1 = high_bumps(w1, p); h2 = high_bumps(w2, p)
                            c = (b2 in h1 and b1 in h2 and h1[b2] == h2[b1])
                            comm[c] += 1
                            if not c and len(ex) < 4:
                                ex.append((word, p, b1, b2, s1, s2))
print(stats)
print("disjoint-support commutation:", comm)
for e in ex:
    print(e)
