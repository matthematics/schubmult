"""Second-order feature search: which features make the join shape deterministic?"""
import sys, itertools
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import rroot, rroots, ltb_down, valid_rc, maxd
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

def shape(word, seq, p, bT, LT):
    BL = {x for x, _ in bumps(LT, seq, p).values()}
    BL2 = set()
    for x in BL:
        BL2 |= {y for y, _ in bumps(x, seq, p).values()}
    cur = bT
    for k in range(4):
        if cur is None:
            return "none"
        if cur == LT:
            return f"L{k}=LT"
        if cur in BL:
            return f"L{k}=b'LT"
        if cur in BL2:
            return f"L{k}=b'b''LT"
        cur, _ = Lchain(cur, p)
    return "none"

def srt(r):
    return (min(r), max(r))

records = []
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
                R = [srt(r) for r in rroots(word)]
                for qb, (bT, Cb) in bumps(word, seq, p).items():
                    if qb == qL:
                        continue
                    a, b = R[qb]
                    wb = perm_of(bT); mb = maxd(wb)
                    LbT, CLb = Lchain(bT, p)
                    # chain of L in beta T vs in T
                    f = {
                        "qb_in_CL": qb in CL,
                        "qL_in_Cb": qL in Cb,
                        "disjoint": not (set(CL) & set(Cb)),
                        "share": "a=m" if a == m else "b=m" if b == m else "b=m+1" if b == m + 1 else "none",
                        "qb>qL": qb > qL,
                        "beta_lowers_maxd": mb < m,
                        "Cb_end_lt_CL_end": min(Cb) < min(CL),
                        "Cb_end_eq_CL_end": min(Cb) == min(CL),
                        "CL_same_in_bT": (CLb == CL) if CLb is not None else None,
                        "CL_bT_meets_Cb": bool(set(CLb) & set(Cb)) if CLb is not None else None,
                        "c_le_p": min(R[min(Cb)]) <= p,  # root of beta's last chain position before bump has coords; the new root (c,a): c = min coordinate of chain end root after? approx
                        "lenCb": min(len(Cb), 3),
                        "lenCL": min(len(CL), 3),
                        "pipes_shared": len({wa[a - 1], wa[b - 1]} & {wa[m - 1], wa[m]}),
                    }
                    records.append((f, shape(word, seq, p, bT, LT), (rows[:p], p, qb)))

print("records", len(records))
keys = list(records[0][0].keys())
def determinism(feat_subset):
    tbl = defaultdict(set)
    for f, s, _ in records:
        tbl[tuple(f[k] for k in feat_subset)].add(s)
    return sum(1 for v in tbl.values() if len(v) > 1), len(tbl)

# greedy forward selection
chosen = []
for step in range(len(keys)):
    best = None
    for k in keys:
        if k in chosen:
            continue
        amb, n = determinism(chosen + [k])
        if best is None or amb < best[0]:
            best = (amb, n, k)
    chosen.append(best[2])
    print(f"+{best[2]}: ambiguous classes {best[0]} of {best[1]}")
    if best[0] == 0:
        break
# print the ambiguous classes for the final selection
tbl = defaultdict(Counter); exs = defaultdict(list)
for f, s, ex in records:
    key = tuple(f[k] for k in chosen)
    tbl[key][s] += 1
    if len(exs[(key, s)]) < 2:
        exs[(key, s)].append(ex)
print("chosen:", chosen)
for key, c in tbl.items():
    if len(c) > 1:
        print(dict(zip(chosen, key)), dict(c))
        for s in c:
            print("     ", s, exs[(key, s)])
