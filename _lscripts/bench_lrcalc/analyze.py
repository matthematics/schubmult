"""Analyze exhaustive S_n timing records from harness.c / bench_schubmult.py.

Usage: python analyze.py n c.bin rs.bin sm.bin
Cross-checks results (nterms + hash) between implementations, then reports timing distributions and
win/loss broken down by structural features of the pair: l(u)+l(v), length of theta(v^{-1}) (v-path depth
for the schubmult kernel on the longer factor), the number of monomials of the shorter factor (lrcalc's
transition expansion size), output size.
"""
import itertools
import sys

import numpy as np

REC = np.dtype([("ui", "<u4"), ("vi", "<u4"), ("t", "<u4"), ("nt", "<u4"), ("h", "<u8")])


def load(path):
    a = np.fromfile(path, dtype=REC)
    order = np.lexsort((a["vi"], a["ui"]))
    return a[order]


def perm_features(n):
    from schubmult import Permutation, Sx, uncode
    from schubmult.symbolic import expand
    from schubmult.utils.schub_lib import compute_vpathdicts

    perms = [Permutation(list(p)) for p in itertools.permutations(range(1, n + 1))]
    inv = np.array([p.inv for p in perms])
    theta_len = np.array([len([t for t in (~p).theta() if t]) for p in perms])
    nmon = np.empty(len(perms), dtype=np.int64)
    nvpath = np.empty(len(perms), dtype=np.int64)  # total v-path states schubmult visits for v
    for i, p in enumerate(perms):
        e = expand(Sx(p).expand())
        nmon[i] = len(e.args) if e.is_Add else 1
        th = [t for t in (~p).theta()]
        while th and th[-1] == 0:
            th.pop()
        if not th:
            nvpath[i] = 0
            continue
        vpd = compute_vpathdicts(th, p * uncode(th))
        nvpath[i] = sum(len(steps) for layer in vpd for steps in layer.values())
    return inv, theta_len, nmon, nvpath


def spearman(x, y):
    rx = np.argsort(np.argsort(x)).astype(np.float64)
    ry = np.argsort(np.argsort(y)).astype(np.float64)
    return np.corrcoef(rx, ry)[0, 1]


def fit_cost(t, feats, names):
    """log t ~ a + sum b_i log(1 + f_i); report coefficients and R^2."""
    X = np.column_stack([np.ones(len(t))] + [np.log1p(f.astype(np.float64)) for f in feats])
    y = np.log(t)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ coef
    r2 = 1 - resid.var() / y.var()
    return " + ".join([f"{coef[0]:.2f}"] + [f"{c:.2f}*log(1+{nm})" for c, nm in zip(coef[1:], names)]), r2


def pct(a, q):
    return np.percentile(a, q)


def fmt_us(ns):
    return f"{ns / 1000:8.1f}"


def table(title, keys, groups, names, arrays):
    print(f"\n{title}")
    hdr = f"{'bucket':>14} {'count':>8} " + " ".join(f"{nm + ' med us':>14}" for nm in names) + "  " + " ".join(f"{nm + ' wins':>9}" for nm in names)
    print(hdr)
    for k in keys:
        m = groups == k
        c = int(m.sum())
        if c == 0:
            continue
        meds = [np.median(a[m]) / 1000 for a in arrays]
        stacked = np.stack([a[m] for a in arrays])
        wins = [(stacked.argmin(axis=0) == i).sum() for i in range(len(arrays))]
        print(f"{str(k):>14} {c:>8} " + " ".join(f"{x:>14.1f}" for x in meds) + "  " + " ".join(f"{w:>9}" for w in wins))


def main():
    n = int(sys.argv[1])
    names = ["lrcalc-C1.2", "lrcalc-rs", "schubmult"] if len(sys.argv) < 6 else ["lrcalc-C1.2", "lrcalc-rs", "sm-native", "sm-python"]
    data = [load(p) for p in sys.argv[2:]]
    N = len(data[0])
    assert all(len(d) == N for d in data), [len(d) for d in data]
    assert all((d["ui"] == data[0]["ui"]).all() and (d["vi"] == data[0]["vi"]).all() for d in data)

    # correctness
    for i in range(1, len(data)):
        bad = ((data[i]["nt"] != data[0]["nt"]) | (data[i]["h"] != data[0]["h"])).sum()
        print(f"{names[i]} vs {names[0]}: {bad} result mismatches out of {N}")

    T = [d["t"].astype(np.float64) for d in data]
    print(f"\nS_{n}: {N} products (u <= v)")
    print(f"{'':>14} {'total s':>9} {'mean us':>9} {'median':>9} {'p90':>9} {'p99':>9} {'p99.9':>9} {'max':>9}")
    for nm, t in zip(names, T):
        print(f"{nm:>14} {t.sum() / 1e9:9.2f} {fmt_us(t.mean())} {fmt_us(pct(t, 50))} {fmt_us(pct(t, 90))} {fmt_us(pct(t, 99))} {fmt_us(pct(t, 99.9))} {fmt_us(t.max())}")

    stacked = np.stack(T)
    winner = stacked.argmin(axis=0)
    print("\nper-pair fastest:", {nm: int((winner == i).sum()) for i, nm in enumerate(names)})
    for i, j in [(k, 1) for k in range(len(T)) if k != 1] + [(1, 0)]:
        r = T[i] / T[j]
        print(f"ratio {names[i]}/{names[j]}: geomean {np.exp(np.log(r).mean()):.3f}, median {np.median(r):.3f}, p10 {pct(r, 10):.3f}, p90 {pct(r, 90):.3f}, min {r.min():.3f}, max {r.max():.1f}")

    inv, theta_len, nmon, nvpath = perm_features(n)
    ui, vi = data[0]["ui"], data[0]["vi"]
    nt = data[0]["nt"].astype(np.int64)
    linv = inv[ui] + inv[vi]
    # lrcalc expands the shorter factor; schubmult's v-path recursion runs over the second argument v
    short = np.where(inv[ui] <= inv[vi], ui, vi)
    nmon_short = nmon[short]
    nmon_min = np.minimum(nmon[ui], nmon[vi])
    th_v = theta_len[vi]
    nvp_v = nvpath[vi]
    nvp_min = np.minimum(nvpath[ui], nvpath[vi])

    table("by l(u)+l(v)", sorted(set(linv)), linv, names, T)
    table("by theta-length of v (schubmult v-path depth)", sorted(set(th_v)), th_v, names, T)
    nb = np.digitize(nmon_short, [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
    labels = {i: f"[{[0,1,2,4,8,16,32,64,128,256,512,1024][i-1]},{[1,2,4,8,16,32,64,128,256,512,1024,99999][i-1]})" for i in range(1, 13)}
    nb_l = np.array([labels.get(b, str(b)) for b in nb])
    table("by #monomials of shorter factor (lrcalc trans() size)", [labels[i] for i in range(1, 13) if labels[i] in set(nb_l)], nb_l, names, T)
    vb = np.digitize(nvp_v, [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
    vb_l = np.array([labels.get(b, str(b)) for b in vb])
    table("by #v-path states of v (schubmult recursion size)", [labels[i] for i in range(1, 13) if labels[i] in set(vb_l)], vb_l, names, T)
    ob = np.digitize(nt, [1, 2, 4, 8, 16, 32, 64, 128, 256])
    ol = {i: f"[{[0,1,2,4,8,16,32,64,128,256][i-1]},{[1,2,4,8,16,32,64,128,256,99999][i-1]})" for i in range(1, 11)}
    ob_l = np.array([ol.get(b, str(b)) for b in ob])
    table("by #output terms", [ol[i] for i in range(1, 11) if ol[i] in set(ob_l)], ob_l, names, T)

    # extremes (use the native kernel when present)
    si = 2
    print(f"\nlargest {names[si]}/lrcalc-rs ratios ({names[si]} slowest relative):")
    r = T[si] / T[1]
    perms = list(itertools.permutations(range(1, n + 1)))
    for k in np.argsort(-r)[:5]:
        print(f"  u={perms[ui[k]]} v={perms[vi[k]]}  sm {T[si][k]/1000:.1f}us rs {T[1][k]/1000:.1f}us C {T[0][k]/1000:.1f}us  nterms={nt[k]} nmon_short={nmon_short[k]} th_v={th_v[k]} nvp_v={nvp_v[k]}")
    print(f"smallest {names[si]}/lrcalc-rs ratios ({names[si]} fastest relative):")
    for k in np.argsort(r)[:5]:
        print(f"  u={perms[ui[k]]} v={perms[vi[k]]}  sm {T[si][k]/1000:.1f}us rs {T[1][k]/1000:.1f}us C {T[0][k]/1000:.1f}us  nterms={nt[k]} nmon_short={nmon_short[k]} th_v={th_v[k]} nvp_v={nvp_v[k]}")

    print("\nSpearman correlation of log(time) with features:")
    feats = [("l(u)+l(v)", linv), ("nmon_short", nmon_short), ("nmon_min", nmon_min), ("theta_len(v)", th_v), ("nvpath(v)", nvp_v), ("nvpath_min", nvp_min), ("nterms", nt)]
    for nm, t in zip(names, T):
        print(f"  {nm:>12}: " + "  ".join(f"{fname}={spearman(np.log(t), f):.2f}" for fname, f in feats))

    print("\nCost model fits  log(t_ns) ~ a + b*log(1+feature)...:")
    for nm, t in zip(names, T):
        for fs in [["nmon_short", "nterms"], ["nvpath(v)", "nterms"], ["nmon_short", "nvpath(v)", "nterms", "l(u)+l(v)"]]:
            fd = dict(feats)
            formula, r2 = fit_cost(t, [fd[f] for f in fs], fs)
            print(f"  {nm:>12} R^2={r2:.3f}: {formula}")

    # hypothetical hybrid: pick by a cheap structural test before computing
    print("\nHybrid dispatch (choose lrcalc-rs when nmon_short <= K, else native schubmult):")
    best_total = np.minimum(T[1], T[si]).sum() / 1e9
    print(f"  oracle (always fastest of the two): {best_total:.2f}s;  lrcalc-rs alone {T[1].sum()/1e9:.2f}s;  {names[si]} alone {T[si].sum()/1e9:.2f}s")
    for K in [1, 2, 4, 8, 12, 16, 24, 32, 64]:
        pick_rs = nmon_short <= K
        tot = np.where(pick_rs, T[1], T[si]).sum() / 1e9
        print(f"  K={K:>3}: total {tot:.2f}s  ({pick_rs.mean()*100:.1f}% to lrcalc-rs)")


if __name__ == "__main__":
    main()
