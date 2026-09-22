"""Root tableau (cell -> (root,row)) under p-high cover bumps: same Q, so compare cellwise.
Hypothesis: along the chain cells c_0 (start),...,c_k, new root at c_i = old root at c_{i-1}, new root at c_0 = (c,a) or reversed."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroots, ltb_down, valid_rc, maxd
from schubmult.utils.perm_utils import is_reduced
from collections import Counter

def eg_Q(word):
    P, Q = [], []
    for t, x in enumerate(word, 1):
        r = 0
        while True:
            if r == len(P):
                P.append([x]); Q.append([t]); break
            row = P[r]
            bigger = [k for k, y in enumerate(row) if y > x]
            if not bigger:
                row.append(x); Q[r].append(t); break
            k = bigger[0]; y = row[k]
            if y == x + 1 and k > 0 and row[k - 1] == x:
                x = y
            else:
                row[k] = x; x = y
            r += 1
    return Q

def cell_of_position(word):
    """Insert reversed word (HY convention): position q of word corresponds to time len-q."""
    Q = eg_Q(tuple(reversed(word)))
    n = len(word)
    cell = {}
    for i, row in enumerate(Q):
        for j, t in enumerate(row):
            cell[n - t] = (i, j)
    return cell

def is_cover(w, a, b):
    arr = list(w) + list(range(len(w) + 1, b + 2))
    if arr[a - 1] < arr[b - 1]:
        return False
    return not any(arr[b - 1] < arr[c - 1] < arr[a - 1] for c in range(a + 1, b))

def srt(r):
    return (min(r), max(r))

st = Counter(); ex = []
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(2, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                cells = cell_of_position(word)
                R = [srt(r) for r in rroots(word)]
                for q in range(len(word)):
                    a, b = R[q]
                    if a <= p or not is_cover(w, a, b):
                        continue
                    nw, chain = ltb_down(word, q, return_chain=True)
                    if nw is None or not valid_rc(nw, seq):
                        continue
                    cells2 = cell_of_position(nw)
                    st[("Q same", cells == cells2)] += 1
                    R2 = [srt(r) for r in rroots(nw)]
                    # non-chain positions: root unchanged?
                    ch = set(chain)
                    st[("nonchain roots unchanged", all(R[k] == R2[k] for k in range(len(word)) if k not in ch))] += 1
                    # chain shift forward: R2[chain[i]] == R[chain[i-1]] ; R2[chain[0]] == (c,a)?
                    fwd = all(R2[chain[i]] == R[chain[i - 1]] for i in range(1, len(chain)))
                    bwd = all(R2[chain[i]] == R[chain[i + 1]] for i in range(len(chain) - 1))
                    st[("chain shift", "fwd" if fwd else "bwd" if bwd else "neither")] += 1
                    if not fwd and not bwd and len(ex) < 4:
                        ex.append((word, p, (a, b), chain, [R[k] for k in chain], [R2[k] for k in chain]))
print(st)
for e in ex:
    print(e)
