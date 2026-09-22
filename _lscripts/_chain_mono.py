"""Chain structure for the last-descent bump on RC graphs: monotone in position? Rows weakly decreasing?"""
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter

def maxd(w):
    return len(w.trimcode)

def ltb_chain(word, index):
    word = [*word]
    chain = []
    while True:
        if word[index] == 1:
            return None, None
        word[index] -= 1
        chain.append(index)
        if is_reduced(word):
            return tuple(word), chain
        index = find_reduced_fail(word, index)
        if index is None:
            return None, None

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        m = maxd(w)
        if m < 2:
            continue
        for Y in RCGraph.all_rc_graphs(w, m):
            rows = [tuple(r) for r in Y]
            if len(rows[-1]) != 0:
                continue
            word, seq = Y.as_reduced_compatible()
            Y1 = Y.little_bump_desc()
            w1, _ = Y1.as_reduced_compatible()
            # find start position: the one whose right-to-left root is (m, m+1)
            starts = []
            for q in range(len(word)):
                new, chain = ltb_chain(word, q)
                if new == tuple(w1):
                    starts.append((q, chain))
            # choose the chain that starts at the position with right root (m,m+1)
            # right-to-left root of position q: s_{a_n}...s_{a_{q+1}} applied to (a_q, a_q+1)
            def rroot(q):
                a, b = word[q], word[q] + 1
                for x in word[q + 1:]:
                    a = a + 1 if a == x else a - 1 if a == x + 1 else a
                    b = b + 1 if b == x else b - 1 if b == x + 1 else b
                return (min(a, b), max(a, b))
            cands = [(q, ch) for q, ch in starts if rroot(q) == (m, m + 1)]
            if len(cands) != 1:
                stats["ambiguous-start"] += 1
                continue
            q, chain = cands[0]
            dec = all(chain[k] > chain[k + 1] for k in range(len(chain) - 1))
            rows_chain = [seq[c] for c in chain]
            rows_dec = all(rows_chain[k] >= rows_chain[k + 1] for k in range(len(chain) - 1))
            stats[("posdec", dec, "rowsweakdec", rows_dec)] += 1
            if not dec and len(bad) < 5:
                bad.append((rows, word, chain, rows_chain))
print(stats)
for b in bad:
    print(b)
