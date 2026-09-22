"""Checks prompted by the external analysis.
(G1) Bumpable positions = Bruhat covers (nearly reduced). Compare with ltb_down degeneracy.
(G2) Entry position of a Z-chain into T is a cover of w_T with a>p.
(G3) No-shift: cover bumps on RC graphs in RC_{<=p} never hit letter 1 / non-unique partner.
(CONF) confluence of cover-only p-high bumps on RC graphs (valid seq).
(LC) Candidate Lemma C: for T, at most one w' in Tset(w_T) with v in P_1(w'), ell(w')-ell(v)=a.
"""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroots, rroot, ltb_down, valid_rc, maxd
from schubmult.utils.perm_utils import is_reduced
from collections import Counter, defaultdict

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])

def is_cover(w, a, b):
    """(a,b) positions, a<b, w(a)>w(b), no c in (a,b) with w(b)<w(c)<w(a)."""
    arr = list(w) + list(range(len(w) + 1, b + 2))
    if arr[a - 1] < arr[b - 1]:
        return False
    return not any(arr[b - 1] < arr[c - 1] < arr[a - 1] for c in range(a + 1, b))

def Dp(X, p):
    h = max(p, maxd(X.perm)); X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

g1 = Counter(); g2 = Counter(); g3 = Counter(); ex3 = []
conf_tbl = {}
memo = {}

def cover_bumps(word, seq, p):
    w = perm_of(word)
    out = set()
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a <= p or not is_cover(w, a, b):
            continue
        nw = ltb_down(word, q)
        if nw is None:
            g3["degenerate (letter1/partner)"] += 1
            if len(ex3) < 3:
                ex3.append((word, seq, q, (a, b)))
            continue
        if not valid_rc(nw, seq):
            g3["invalid rc"] += 1
            continue
        g3["ok"] += 1
        out.add(nw)
    return out

def NF(word, seq, p):
    key = (word, seq, p)
    if key in memo:
        return memo[key]
    nb = cover_bumps(word, seq, p)
    res = frozenset([word]) if not nb else frozenset().union(*(NF(x, seq, p) for x in nb))
    memo[key] = res
    return res

lc = Counter(); exlc = []

def Tn(z):
    """T_n(z) via transition: z' with z ->(n) z', same length. Use zero_out on all RC graphs? Use LS transition
    on permutations: z' = z t_{ns} t_{in} covers with maxd<n after full iteration. Compute via RC graphs instead."""
    raise NotImplementedError

for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        # G1: for every reduced word, position q: cover iff ltb_down defined (word level, ignoring letter 1)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                for q in range(len(word)):
                    a, b = rroot(word, q)
                    cov = is_cover(w, a, b)
                    nw = ltb_down(word, q)
                    g1[("cover", cov, "ltb defined", nw is not None, "letter1" if word[q] == 1 else "")] += 1
                # CONF
                nf = NF(word, seq, p)
                conf_tbl[len(nf)] = conf_tbl.get(len(nf), 0) + 1
        # G2: Z-chain entry into T
        for m in (md,):
            if m < 3:
                continue
            for Y in RCGraph.all_rc_graphs(w, m):
                rows = [tuple(r) for r in Y]
                if len(rows[-1]) != 0:
                    continue
                Y1 = Y.little_bump_desc()
                for p in range(2, m):
                    u = tuple(x for r in rows[:p] for x in r)
                    u1 = tuple(x for r in [tuple(r) for r in Y1][:p] for x in r)
                    if u == u1:
                        continue
                    changed = [k for k in range(len(u)) if u[k] != u1[k]]
                    q = max(changed)
                    a, b = rroot(u, q)
                    g2[("entry cover", is_cover(perm_of(u), a, b), "a>p", a > p)] += 1
print("G1:", g1)
print("G2:", g2)
print("G3:", g3)
for e in ex3:
    print("   ", e)
print("CONF normal form counts:", conf_tbl)
