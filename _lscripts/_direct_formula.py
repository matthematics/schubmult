"""Hypothesis H: D^p(T) == remove all p-high crossings (cover order along last descents), then pieri_insert(p, rows)."""
import sys
sys.path.insert(0, "/home/matthematics/schubmult/_lscripts")
from itertools import permutations
from schubmult import Permutation, RCGraph
from _lb_common import maxd
from collections import Counter

def Dp(X, p):
    h = max(p, maxd(X.perm))
    X = X.resize(h)
    while len(X) > p:
        X = X.zero_out_last_row()
    return X

def flat(T, p):
    """Remove crossings at root (m,m+1), m = current maxd, until maxd <= p. Return (Tflat, rows, descs)."""
    interim = T
    rows = []; descs = []
    while maxd(interim.perm) > p:
        m = maxd(interim.perm)
        descs.append(m)
        interim, row = interim.exchange_property(m, return_row=True)
        rows.append(row)
    return interim, rows, descs

def H_variants(T, p):
    Tf, rows, descs = flat(T, p)
    out = {}
    # V1: plain insertion at height p (no placeholder)
    try:
        X = Tf.resize(p).pieri_insert(p, rows) if rows else Tf.resize(p)
        out["V1"] = X.resize(p) if len(X) >= p else None
    except Exception as e:
        out["V1"] = ("err", str(e)[:40])
    # V2: placeholder row p+1 with letters descs (like zero_out_last_row), insert at descent p, then drop row p+1
    try:
        base = Tf.resize(p + 1)
        interim2 = RCGraph([*base[:p], tuple(sorted(descs, reverse=True))])
        X = interim2.pieri_insert(p, rows) if rows else interim2
        out["V2"] = X.rowrange(0, p)
    except Exception as e:
        out["V2"] = ("err", str(e)[:40])
    # V3: placeholder letters p+1..p+len(rows)
    try:
        base = Tf.resize(p + 1)
        ph = tuple(range(p + len(rows), p, -1))
        interim2 = RCGraph([*base[:p], ph])
        X = interim2.pieri_insert(p, rows) if rows else interim2
        out["V3"] = X.rowrange(0, p)
    except Exception as e:
        out["V3"] = ("err", str(e)[:40])
    return out

stats = Counter()
bad = {}
for N in range(2, 8):
    for arr in permutations(range(1, N + 1)):
        w = Permutation(list(arr))
        if w.inv > 7 or w.inv == 0:
            continue
        md = maxd(w)
        for p in range(1, md):
            if md - p < 2:
                continue  # single level is the algorithm itself; test multi-level
            for T in RCGraph.all_rc_graphs(w, md):
                rows = [tuple(r) for r in T]
                if any(len(rows[i]) for i in range(p, md)):
                    continue
                D = Dp(T, p)
                for k, v in H_variants(T, p).items():
                    ok = (v == D)
                    stats[(k, ok)] += 1
                    if not ok and k not in bad:
                        bad[k] = (rows, p, [tuple(r) for r in D], v if isinstance(v, tuple) else [tuple(r) for r in v] if v is not None else None)
print(stats)
for k, b in bad.items():
    print(k, b)
