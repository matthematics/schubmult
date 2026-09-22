"""Correct emulation of zeta on Y = T u row_{p+1}(c) (padded to height m = maxd(w)):
iterate the Little bump at the root alpha_{m'} of the *current* permutation (m' = its maxd)
until maxd < m.  Check (1) emulation == zero_out_last_row on the top p rows; (2) every bump
stays inside the word of T and equals the bump of word(T) alone at the same position (Fact A);
(3) record the sequence of a-frame roots sigma_c(alpha_{m'}) and compare with clip's own chain.
"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from collections import Counter
from schubmult import Permutation, RCGraph
from schubmult.utils.perm_utils import is_reduced
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

def zeta_emulate(word, la, m):
    """Iterate bump at alpha_{maxd} on the full word until maxd < m. Returns (final word, info)."""
    steps = []  # (m', start position, chain, a-frame root)
    while True:
        mp = maxd(perm_of(word))
        if mp < m:
            return word, steps
        k = pos_of_root(word, (mp, mp + 1))
        if k is None:
            return None, steps
        aroot = rroot(word[:la], k) if k < la else None
        new, chain = ltb_down(word, k, return_chain=True)
        if new is None:
            return None, steps
        steps.append((mp, k, tuple(chain), aroot))
        word = new

def clip_chain(word, p):
    chain = [word]; roots = []
    while True:
        d = maxd(perm_of(word))
        if d <= p:
            return chain, roots
        k = pos_of_root(word, (d, d + 1))
        if k is None:
            return None, None
        new = ltb_down(word, k)
        if new is None:
            return None, None
        roots.append((d, d + 1)); word = new; chain.append(word)

stats = Counter(); examples = {}
MAXLEN = 8
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
                for c in range(1, MAXLEN - w1.inv + 1):
                    sigma = Permutation([])
                    for j in range(p + c, p, -1):
                        sigma = sigma * Permutation.ref_product(j)
                    w = w1 * sigma
                    if w.inv != w1.inv + c:
                        continue
                    m = maxd(w)
                    if m < p + 2:
                        continue
                    case = "A" if c <= m - p - 2 else ("B" if c == m - p - 1 else "C")
                    Yrows = Trows + [tuple(range(p + c, p, -1))] + [()] * (m - p - 1)
                    Y = RCGraph(Yrows)
                    Z = Y.zero_out_last_row(); Zrows = [tuple(r) for r in Z]
                    topZ = Zrows[:p]
                    aw, aseq = word_seq(Trows); la = len(aw)
                    fw, fseq = word_seq(Yrows)
                    final, steps = zeta_emulate(fw, la, m)
                    if final is None:
                        stats[(case, "emulation degenerate")] += 1; continue
                    emu_top = rows_from(final[:la], aseq, p)
                    stats[(case, "emulation == zeta", emu_top == topZ)] += 1
                    if emu_top != topZ:
                        examples.setdefault(("emu mismatch", case), (Yrows, topZ, emu_top, steps)); continue
                    # Fact A: all bump chains inside a; a-alone bumps agree
                    inside = all(k < la and all(i < la for i in ch) for (_, k, ch, _) in steps)
                    stats[(case, "all bumps inside a", inside)] += 1
                    if inside:
                        aw2 = aw; ok = True
                        for (_, k, ch, _) in steps:
                            new = ltb_down(aw2, k)
                            if new is None:
                                ok = False; break
                            aw2 = new
                        stats[(case, "a-alone bumps == restriction", ok and aw2 == final[:la])] += 1
                    # roots used, in a-frame, and maxd excursion
                    aroots = tuple(s[3] for s in steps)
                    mprimes = tuple(s[0] for s in steps)
                    stats[(case, "maxd rises above m", any(mp > m for mp in mprimes))] += 1
                    # compare with clip chain on a
                    cchain, croots = clip_chain(aw, p)
                    zchain = [aw]; aw2 = aw
                    for (_, k, ch, _) in steps:
                        aw2 = ltb_down(aw2, k); zchain.append(aw2)
                    prefix = cchain is not None and cchain[:len(zchain)] == zchain
                    stats[(case, "zeta chain prefix of clip chain", prefix)] += 1
                    if not prefix:
                        key = (case, "nonprefix", aroots[:2], tuple(croots[:2]) if croots else None)
                        stats[key] += 1
                        examples.setdefault(key, (Trows, p, c, m, aroots, croots, zchain, cchain))
                    cz = RCGraph(Zrows).__class__  # placeholder
for k in sorted(stats, key=str):
    print(k, stats[k])
print()
for k, e in examples.items():
    print(k, e)
