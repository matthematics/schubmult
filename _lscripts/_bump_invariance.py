"""Test (L): T with rows > p empty (height h = max(p, maxd)), q a position whose left-to-right root (a,b) has a > p,
and ltb_q(T) (downward word bump, same compatible sequence, valid RC letters >= row) => D^p(ltb_q T) == D^p(T).
Compare with positions having a <= p."""
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

def Dp(X, p):
    rows = [tuple(r) for r in X][:p]
    top = RCGraph(rows)
    h = max(p, maxd(top.perm))
    top = top.resize(h)
    while len(top) > p:
        top = top.zero_out_last_row()
    return top

def valid_rc(word, seq):
    # letters >= row, strictly decreasing within a row
    for k in range(len(word)):
        if word[k] < seq[k]:
            return False
        if k + 1 < len(word) and seq[k] == seq[k + 1] and not word[k] > word[k + 1]:
            return False
    return True

stats = Counter()
bad = []
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md + 1):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible()
                DT = Dp(T, p)
                for q in range(len(word)):
                    a, b = T.left_to_right_inversion(q)
                    nw = ltb_down(word, q)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    T2 = RCGraph.from_reduced_compatible(nw, seq)
                    ok = Dp(T2, p) == DT
                    stats[("a>p" if a > p else "a<=p", ok)] += 1
                    if a > p and not ok and len(bad) < 5:
                        bad.append((rows, p, q, (a, b), nw))
print(stats)
for b in bad:
    print(b)
