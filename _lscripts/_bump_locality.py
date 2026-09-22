"""Locality of p-high bumps: (i) non-chain positions keep their roots; (ii) chain roots are (a,b0),(a,b1),..,(c,a)
with b0>b1>...>a>c; (iii) low-low subword invariant; (iv) after bump, chain positions carry roots
(a,b1),...,(a,b_{k-1}),(c,a)."""
import random, sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation
from _lb_common import rroots, reduced_words, maxd
from schubmult.utils.perm_utils import is_reduced
from collections import Counter

def chain(word, index):
    word = [*word]; ch = []
    while True:
        if word[index] == 1:
            return None
        ch.append(index)
        word[index] -= 1
        if is_reduced(word):
            return tuple(word), ch
        rts = rroots(word)
        target = set(rts[index])
        partners = [i for i in range(len(word)) if i != index and set(rts[i]) == target]
        if len(partners) != 1:
            return None
        index = partners[0]

def srt(r):
    return (min(r), max(r))

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
            rts = [srt(r) for r in rroots(word)]
            for p in range(1, maxd(w)):
                for q in range(len(word)):
                    a, b = rts[q]
                    if a <= p:
                        continue
                    res = chain(word, q)
                    if res is None:
                        continue
                    nw, ch = res
                    nrts = [srt(r) for r in rroots(nw)]
                    chs = set(ch)
                    stats[("(i) nonchain roots fixed", all(rts[k] == nrts[k] for k in range(len(word)) if k not in chs))] += 1
                    roots_before = [rts[k] for k in ch]
                    bs = [r[1] for r in roots_before[:-1]]
                    last = roots_before[-1]
                    ok2 = all(r[0] == a for r in roots_before[:-1]) and all(bs[i] > bs[i + 1] for i in range(len(bs) - 1)) and (len(ch) == 1 or (last[1] == a and last[0] < a) or (last[0] == a))
                    stats[("(ii) chain roots (a,b_i) decreasing", ok2)] += 1
                    # last root: (c,a) with c<a, or if chain length 1 the single root is (a,b) itself
                    roots_after = [nrts[k] for k in ch]
                    shift_ok = roots_after[:-1] == roots_before[1:] if len(ch) > 1 else True
                    stats[("(iv) roots shift along chain", shift_ok)] += 1
                    ll = tuple(word[k] for k in range(len(word)) if rts[k][1] <= p)
                    nll = tuple(nw[k] for k in range(len(word)) if nrts[k][1] <= p)
                    stats[("(iii) low-low subword fixed", ll == nll)] += 1
                    if not shift_ok and len(ex) < 4:
                        ex.append((word, p, q, ch, roots_before, roots_after))
print(stats)
for e in ex:
    print(e)
