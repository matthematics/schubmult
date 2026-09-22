"""Does a small set of chain-interaction features determine the join shape?
Features: qb in C(L)? qL in C(beta)? chains disjoint? shared pipe? a==m? position order of starts."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, ltb_down, valid_rc, maxd
from collections import Counter, defaultdict

def perm_of(word):
    return Permutation.ref_product(*word) if word else Permutation([])

def arr(w, n=14):
    a = list(w); a += list(range(len(a) + 1, n + 1)); return a[:n]

def is_cover(w, a, b):
    ar = arr(w)
    if ar[a - 1] < ar[b - 1]:
        return False
    return not any(ar[b - 1] < ar[c - 1] < ar[a - 1] for c in range(a + 1, b))

def Lchain(word, p):
    w = perm_of(word); m = maxd(w)
    if m <= p:
        return None, None
    q = next(q for q in range(len(word)) if rroot(word, q) == (m, m + 1))
    return ltb_down(word, q, return_chain=True)

def bumps(word, seq, p):
    w = perm_of(word); out = {}
    for q in range(len(word)):
        a, b = rroot(word, q)
        if a > p and is_cover(w, a, b):
            res = ltb_down(word, q, return_chain=True)
            if res[0] is not None and valid_rc(res[0], seq):
                out[q] = res
    return out

def shape(word, seq, p, bT):
    LT, _ = Lchain(word, p)
    BL = {x for x, _ in bumps(LT, seq, p).values()}
    BL2 = set()
    for x in BL:
        BL2 |= {y for y, _ in bumps(x, seq, p).values()}
    cur = bT
    for k in range(4):
        if cur is None:
            return "none"
        if cur == LT:
            return f"L^{k}(bT)=LT"
        if cur in BL:
            return f"L^{k}(bT)=b'(LT)"
        if cur in BL2:
            return f"L^{k}(bT)=b'b''(LT)"
        cur, _ = Lchain(cur, p)
    return "none"

tbl = defaultdict(Counter)
for N in range(2, 8):
    for a_ in permutations(range(1, N + 1)):
        w = Permutation(list(a_))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w); m = md
        wa = arr(w)
        for p in range(1, md):
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                word, seq = T.as_reduced_compatible(); word = tuple(word); seq = tuple(seq)
                LT, CL = Lchain(word, p)
                if LT is None:
                    continue
                qL = CL[0]
                for qb, (bT, Cb) in bumps(word, seq, p).items():
                    if qb == qL:
                        continue
                    a, b = rroot(word, qb)
                    feats = (
                        "qb in CL" if qb in CL else "qb notin CL",
                        "qL in Cb" if qL in Cb else "qL notin Cb",
                        "disjoint" if not (set(CL) & set(Cb)) else "overlap",
                        "a==m" if a == m else ("b==m" if b == m else ("b==m+1" if b == m + 1 else "no shared pos")),
                        "qb>qL" if qb > qL else "qb<qL",
                    )
                    tbl[feats][shape(word, seq, p, bT)] += 1
for f, c in sorted(tbl.items(), key=lambda x: -sum(x[1].values())):
    print(sum(c.values()), f, dict(c))
