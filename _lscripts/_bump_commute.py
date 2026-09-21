"""Test commutation: for T (rows>p empty, maxd = h > p) and a valid bump q with root a>p:
lmap(ltb_q T) == lmap(T) or == ltb_{q'}(lmap T) for some q' with root a'>p."""
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter

def maxd(w):
    return len(w.trimcode)

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

def big_bumps(T, p):
    word, seq = T.as_reduced_compatible()
    out = {}
    for q in range(len(word)):
        a, b = T.left_to_right_inversion(q)
        if a <= p:
            continue
        nw = ltb_down(word, q)
        if nw is None or not valid_rc(nw, seq):
            continue
        out[q] = RCGraph.from_reduced_compatible(nw, seq)
    return out

def lmap(T, h):
    T = T.resize(h)
    return T.little_bump_desc()

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        h = maxd(w)
        for p in range(2, h):
            for T in RCGraph.all_rc_graphs(w, h):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, h)):
                    continue
                B = lmap(T, h)
                Bbumps = {tuple(tuple(r) for r in X.resize(h)) for X in big_bumps(B, p).values()}
                Bkey = tuple(tuple(r) for r in B.resize(h))
                for q, T2 in big_bumps(T, p).items():
                    if maxd(T2.perm) < h:
                        stats["bumped-below-h"] += 1
                        continue
                    A = lmap(T2, h)
                    Akey = tuple(tuple(r) for r in A.resize(h))
                    if Akey == Bkey:
                        stats["equal"] += 1
                    elif Akey in Bbumps:
                        stats["commute"] += 1
                    else:
                        stats["FAIL"] += 1
                        if len(bad) < 5:
                            bad.append((rows, p, q, [tuple(r) for r in T2], Akey, Bkey))
print(stats)
for b in bad:
    print(b)
