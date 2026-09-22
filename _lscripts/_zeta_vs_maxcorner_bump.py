"""Is zero_out_last_row equal to iterated Little bumps at the MAXIMAL-CORNER root (r,s)
(Little's original transition bump / inverse of Gao-Huang's m_{i,r} on pipe dreams),
iterated while maxd == n?  Also check each step is a transition step w -> w t_{rs} t_{ir}."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from collections import Counter
from schubmult import Permutation, RCGraph
from _lb_common import ltb_down, rroot, maxd

def perm_of(word):
    w = Permutation([])
    for a in word:
        w = w * Permutation.ref_product(a)
    return w

def word_seq(rows):
    word, seq = [], []
    for i, r in enumerate(rows, 1):
        for a in sorted(r, reverse=True):
            word.append(a); seq.append(i)
    return tuple(word), tuple(seq)

def rows_from(word, seq, height):
    rows = [[] for _ in range(height)]
    for a, i in zip(word, seq):
        rows[i - 1].append(a)
    return [tuple(sorted(r, reverse=True)) for r in rows]

def max_corner(w):
    r = maxd(w)
    s = max(j for j in range(r + 1, len(w) + 2) if (w[j - 1] if j - 1 < len(w) else j) < w[r - 1])
    return r, s

def is_transition_step(w, w2, r, s):
    x = w * Permutation.ref_product  # placeholder
    return None

stats = Counter(); ex = []
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5
for N in range(2, NMAX + 1):
    for a in permutations(range(1, N + 1)):
        w = Permutation(list(a))
        if w.inv == 0:
            continue
        d = maxd(w)
        for n in range(max(d, 2), N + 1):
            for R in RCGraph.all_rc_graphs(w, n):
                if len(R[-1]) != 0:
                    continue
                Z = R.zero_out_last_row(); Zrows = [tuple(r) for r in Z]
                word, seq = word_seq([tuple(r) for r in R])
                cur = word; ok = True; steps = 0
                while maxd(perm_of(cur)) >= n:
                    wc = perm_of(cur)
                    r, s = max_corner(wc)
                    k = [i for i in range(len(cur)) if rroot(cur, i) == (r, s)]
                    if len(k) != 1:
                        ok = False; break
                    new = ltb_down(cur, k[0])
                    if new is None:
                        ok = False; break
                    wn = perm_of(new)
                    # transition step check: wn = wc t_{rs} t_{ir} for some i<r, length preserved
                    wt = wc.swap(r - 1, s - 1)
                    good = any(wt.swap(i - 1, r - 1) == wn for i in range(1, r)) and wn.inv == wc.inv
                    stats[("step is transition", good)] += 1
                    stats[("maxd after step <= n", maxd(wn) <= n)] += 1
                    cur = new; steps += 1
                if not ok:
                    stats["degenerate"] += 1; continue
                emu = rows_from(cur, seq, n - 1)
                eq = emu == Zrows
                stats[("maxcorner emulation == zeta", eq)] += 1
                if not eq and len(ex) < 5:
                    ex.append(([tuple(r) for r in R], Zrows, emu))
for k, v in sorted(stats.items(), key=str):
    print(v, k)
for e in ex:
    print("  ", e)
