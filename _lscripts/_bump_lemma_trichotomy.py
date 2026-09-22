"""Bump lemma / trichotomy check for clpzero in the padded single-row setting.

Y = T u row_{p+1}(c), padded to height m = maxd(w), w = w_1 sigma_c, sigma_c = s_{p+c}...s_{p+1}.
Claim: zeta acts on T by iterating the Little bump of word(T) alone at the root sigma_c(alpha_m),
stopping when maxd(w_1^{(k)} sigma_c) < m.  Cases: A: c<=m-p-2 (root alpha_m), B: c=m-p-1 (root (m-1,m+1)),
C: c>=m-p (root alpha_{m-1}).  Also compare zeta's bump sequence with clip's own bump sequence on T.
"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from collections import Counter
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd

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

from schubmult.utils.perm_utils import is_reduced
not_cover = Counter()

def bump_at_root(word, seq, root):
    """Little bump of word at the unique position with right root == root; None if absent/degenerate.

    Records in not_cover whether deleting the bump letter fails to leave a reduced word."""
    pos = [k for k in range(len(word)) if rroot(word, k) == root]
    if len(pos) != 1:
        return None
    k = pos[0]
    not_cover[is_reduced(list(word[:k]) + list(word[k + 1:]))] += 1
    new, chain = ltb_down(word, k, return_chain=True)
    if new is None or not valid_rc(new, seq):
        return None
    return new

def clip_direct(rows, p):
    T = RCGraph([tuple(r) for r in rows[:p]] + [()] * 0)
    h = max(p, maxd(T.perm)); T = T.resize(h)
    while len(T) > p:
        T = T.zero_out_last_row()
    return [tuple(r) for r in T]

def clip_chain(word, seq, p):
    """clip's own bump sequence on word: bump at alpha_{maxd} while maxd > p. Returns list of words."""
    chain = [word]
    while True:
        w1 = perm_of(word); d = maxd(w1)
        if d <= p:
            return chain
        new = bump_at_root(word, seq, (d, d + 1))
        if new is None:
            return None
        word = new; chain.append(word)

stats = Counter(); bad = []; caseB = []; caseC = []
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
                    Yrows = Trows + [tuple(range(p + c, p, -1))] + [()] * (m - p - 1)
                    Y = RCGraph(Yrows)
                    assert Y.perm == w and len(Y) == m
                    Z = Y.zero_out_last_row()
                    Zrows = [tuple(r) for r in Z]
                    assert Zrows[p] == tuple(range(p + c, p, -1)), ("row p+1 changed", Yrows, Zrows)
                    topZ = Zrows[:p]
                    # predicted: iterate bump at sigma_c(alpha_m) on word(T)
                    beta = tuple(sorted((sigma[m - 1], sigma[m])))
                    if c <= m - p - 2:
                        case = "A"; assert beta == (m, m + 1)
                    elif c == m - p - 1:
                        case = "B"; assert beta == (m - 1, m + 1)
                    else:
                        case = "C"; assert beta == (m - 1, m)
                    word, seq = word_seq(Trows)
                    zchain = [word]; ok = True
                    while maxd(perm_of(word) * sigma) >= m:
                        new = bump_at_root(word, seq, beta)
                        if new is None:
                            ok = False; break
                        word = new; zchain.append(word)
                    pred = rows_from(word, seq, p) if ok else None
                    match = ok and pred == topZ
                    stats[(case, "bump lemma", match)] += 1
                    if not match and len(bad) < 5:
                        bad.append((case, Yrows, topZ, pred))
                    # clpzero
                    cz = clip_direct(Zrows, p) == clip_direct(Yrows, p)
                    stats[(case, "clpzero", cz)] += 1
                    # compare with clip's chain on T
                    cchain = clip_chain(*word_seq(Trows), p)
                    if cchain is None:
                        stats[(case, "clip chain degenerate")] += 1
                        continue
                    prefix = len(zchain) <= len(cchain) and cchain[:len(zchain)] == zchain
                    stats[(case, "zeta chain is prefix of clip chain", prefix)] += 1
                    if case == "B":
                        # is bump at (m-1,m+1) == bump at alpha_m then alpha_{m-1}?
                        w0, s0 = word_seq(Trows)
                        b1 = bump_at_root(w0, s0, (m - 1, m + 1))
                        x = bump_at_root(w0, s0, (m, m + 1))
                        b2 = bump_at_root(x, s0, (m - 1, m)) if x is not None else None
                        x2 = bump_at_root(w0, s0, (m - 1, m))
                        b3 = bump_at_root(x2, s0, (m, m + 1)) if x2 is not None else None
                        stats[("B", "t_{m-1,m+1} == a_m then a_{m-1}", b1 is not None and b1 == b2)] += 1
                        stats[("B", "t_{m-1,m+1} == a_{m-1} then a_m", b1 is not None and b1 == b3)] += 1
                        stats[("B", "maxd(w1)", maxd(w1) - m)] += 1
                        if not prefix and len(caseB) < 4:
                            caseB.append((Trows, p, c, m, zchain, cchain))
                    if case == "C":
                        stats[("C", "maxd(w1)==m-1", maxd(w1) == m - 1, "prefix", prefix)] += 1
                        if not prefix and len(caseC) < 4:
                            caseC.append((Trows, p, c, m, zchain, cchain))
for k in sorted(stats, key=str):
    print(k, stats[k])
print("bump position deletion leaves reduced word (True/False counts):", dict(not_cover))
print("bump lemma failures:", bad)
print("case B non-prefix examples:")
for e in caseB: print("  ", e)
print("case C non-prefix examples:")
for e in caseC: print("  ", e)
