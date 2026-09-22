"""clpzero with an arbitrary (not left-justified) row p+1.

Y = T u rho, rho = decreasing letters >= p+1 in row p+1, padded to height m = maxd(w), w = w1*u,
u = perm(rho), length-additive, m >= p+2.  Emulate zeta by bumps at alpha_{maxd(current)} until maxd < m.
Record: clpzero; whether row p+1 changes; whether bumps start/enter rho; where the zeta-chain on the
top part and clip's chain on T merge (number of clip steps on each side needed to meet).
"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations, combinations
from collections import Counter
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, maxd

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

def perm_of(word):
    w = Permutation([])
    for a in word:
        w = w * Permutation.ref_product(a)
    return w

def pos_of_root(word, root):
    pos = [k for k in range(len(word)) if rroot(word, k) == root]
    return pos[0] if len(pos) == 1 else None

def zeta_emulate(word, m):
    steps = []
    while True:
        mp = maxd(perm_of(word))
        if mp < m:
            return word, steps
        k = pos_of_root(word, (mp, mp + 1))
        if k is None:
            return None, steps
        new, chain = ltb_down(word, k, return_chain=True)
        if new is None:
            return None, steps
        steps.append((mp, k, tuple(chain)))
        word = new

def clip_chain(word, p, limit=30):
    chain = [word]
    while len(chain) < limit:
        d = maxd(perm_of(word))
        if d <= p:
            return chain
        k = pos_of_root(word, (d, d + 1))
        if k is None:
            return None
        word = ltb_down(word, k)
        if word is None:
            return None
        chain.append(word)
    return None

def clip_direct(rows, p):
    T = RCGraph([tuple(r) for r in rows[:p]])
    h = max(p, maxd(T.perm)); T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    return [tuple(r) for r in T]

stats = Counter(); examples = {}
MAXLEN = 8; MAXLETTER_OFF = 6
for N in range(2, 8):
    for a in permutations(range(1, N + 1)):
        w1 = Permutation(list(a))
        if w1.inv == 0 or w1.inv > 6:
            continue
        d1 = maxd(w1)
        for p in range(2, 5):
            n = max(p, d1)
            for T in RCGraph.all_rc_graphs(w1, n):
                Trows = [tuple(r) for r in T]
                if any(len(r) for r in Trows[p:]):
                    continue
                Trows = Trows[:p]
                letters = range(p + 1, p + 1 + MAXLETTER_OFF)
                for c in range(1, MAXLEN - w1.inv + 1):
                    for rho in combinations(sorted(letters, reverse=True), c):
                        u = perm_of(rho)
                        if u.inv != c:
                            continue
                        w = w1 * u
                        if w.inv != w1.inv + c:
                            continue
                        m = maxd(w)
                        if m < p + 2:
                            continue
                        leftjust = rho == tuple(range(p + c, p, -1))
                        Yrows = Trows + [tuple(rho)] + [()] * (m - p - 1)
                        Y = RCGraph(Yrows)
                        Z = Y.zero_out_last_row(); Zrows = [tuple(r) for r in Z]
                        topZ = Zrows[:p]; rowZ = Zrows[p]
                        alpha_m_in_rho = maxd(u) >= m and (m in [d + 1 for d in u.descents()])
                        key = ("leftjust", leftjust, "alpha_m in N(u)", alpha_m_in_rho)
                        cz = clip_direct(Zrows, p) == clip_direct(Yrows, p)
                        stats[key + ("clpzero", cz)] += 1
                        stats[key + ("row p+1 changed", rowZ != tuple(rho))] += 1
                        fw, fseq = word_seq(Yrows); la = sum(len(r) for r in Trows)
                        final, steps = zeta_emulate(fw, m)
                        if final is None:
                            stats[key + ("emulation degenerate",)] += 1; continue
                        emu_top = rows_from(final[:la], fseq[:la], p)
                        stats[key + ("emulation == zeta", emu_top == topZ)] += 1
                        touched_rho = any(k >= la or any(i >= la for i in ch) for (_, k, ch) in steps)
                        started_in_rho = any(k >= la for (_, k, ch) in steps)
                        stats[key + ("bumps touch rho", touched_rho, "start in rho", started_in_rho)] += 1
                        # merge distance of clip chains from T and from topZ
                        aw, aseq = word_seq(Trows)
                        zw, _ = word_seq(topZ)
                        c1 = clip_chain(aw, p); c2 = clip_chain(zw, p)
                        if c1 is None or c2 is None:
                            stats[key + ("clip chain degenerate",)] += 1; continue
                        merge = None
                        for i, x in enumerate(c1):
                            if x in c2:
                                merge = (i, c2.index(x)); break
                        stats[key + ("merge (steps from T, steps from zetaTop)", merge)] += 1
                        if merge is not None and merge[0] >= 2 and ("ex", key, merge) not in examples:
                            examples[("ex", key, merge)] = (Trows, rho, p, m, topZ, c1, c2)
for k in sorted(stats, key=str):
    print(k, stats[k])
print()
for k, e in list(examples.items())[:12]:
    print(k, e)
