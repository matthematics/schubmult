"""RC-level confluence test: p-high downward Little bumps that keep the compatible sequence valid.
Normal forms unique? (Expect normal forms = graphs with maxd <= p.)"""
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter

def maxd(w):
    return len(w.trimcode)

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

def valid_rc(word, seq):
    for k in range(len(word)):
        if word[k] < seq[k]:
            return False
        if k + 1 < len(word) and seq[k] == seq[k + 1] and not word[k] > word[k + 1]:
            return False
    return True

def high_bumps(word, seq, p):
    out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p:
            nw = ltb_down(word, q)
            if nw is not None and valid_rc(nw, seq):
                out.add(nw)
    return out

memo = {}
def normal_forms(word, seq, p):
    key = (word, seq, p)
    if key in memo:
        return memo[key]
    nb = high_bumps(word, seq, p)
    if not nb:
        res = frozenset([word])
    else:
        res = frozenset().union(*(normal_forms(x, seq, p) for x in nb))
    memo[key] = res
    return res

stats = Counter()
nfstats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible()
                nf = normal_forms(word, seq, p)
                stats[len(nf)] += 1
                for x in nf:
                    nfstats[maxd(Permutation.ref_product(*x)) <= p] += 1
                if len(nf) > 1 and len(bad) < 5:
                    bad.append((rows, p, sorted(nf)))
print("normal form counts:", stats)
print("normal forms have maxd<=p:", nfstats)
for b in bad:
    print(b)
