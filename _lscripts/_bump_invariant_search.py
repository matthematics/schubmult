"""Find the invariant: is D^p(T) determined by (Q-tableau of word, weight, positions of values 1..p in perm)?
Also: do bumps with a > p preserve positions of values 1..p?"""
from itertools import permutations
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import find_reduced_fail, is_reduced
from collections import Counter, defaultdict

def maxd(w):
    return len(w.trimcode)

def eg_Q(word):
    """Edelman-Greene insertion recording tableau (row insertion with EG rule), as tuple of rows."""
    P = []
    Q = []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]
            y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                # EG: x already present, bump x+1 without changing the row
                x = y
            else:
                row[k] = x
                x = y
            r += 1
    return tuple(tuple(r) for r in Q)

def Dp(X, p):
    rows = [tuple(r) for r in X][:p]
    top = RCGraph(rows)
    h = max(p, maxd(top.perm))
    top = top.resize(h)
    while len(top) > p:
        top = top.zero_out_last_row()
    return top

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

def inv_positions(w, p):
    arr = list(w)
    arr = arr + list(range(len(arr) + 1, max(len(arr), p) + 2))
    return tuple(arr.index(v) + 1 for v in range(1, p + 1))

tbl = defaultdict(set)
tbl2 = defaultdict(set)
pres = Counter()
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
                wt = tuple(len(r) for r in rows[:p])
                DT = Dp(T, p)
                key = (p, eg_Q(word), wt, inv_positions(w, p))
                tbl[key].add(DT)
                key2 = (p, eg_Q(word), wt, tuple(sorted(inv_positions(w, p))))
                tbl2[key2].add(DT)
                for q in range(len(word)):
                    a, b = T.left_to_right_inversion(q)
                    if a <= p:
                        continue
                    nw = ltb_down(word, q)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    w2 = RCGraph.from_reduced_compatible(nw, seq).perm
                    pres[inv_positions(w, p) == inv_positions(w2, p)] += 1
print("bumps a>p preserve positions of values<=p:", pres)
amb = {k: v for k, v in tbl.items() if len(v) > 1}
print("keys", len(tbl), "ambiguous (Q,wt,positions of 1..p):", len(amb))
for k, v in list(amb.items())[:3]:
    print(k, [[tuple(r) for r in x] for x in v])
amb2 = {k: v for k, v in tbl2.items() if len(v) > 1}
print("ambiguous (Q,wt,set of positions):", len(amb2))
