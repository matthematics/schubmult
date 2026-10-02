"""Generate pair families for the scaling experiment and summarize scaling_pairs output.

  python scaling.py gen  > pairs.txt        (random dense, Grassmannian, many-descent families, n = 6..12)
  python scaling.py summarize < results.txt
"""
import itertools
import random
import sys

import numpy as np


def grassmannian(n, k):
    """1-line form with a single descent at k and shape as full as possible: odd/even interleave."""
    return list(range(1, n + 1, 2)) + list(range(2, n + 1, 2)) if k is None else sorted(range(1, n + 1))[:0]


def families(n, rng, count):
    out = []
    perms_sample = lambda: rng.sample(range(1, n + 1), n)
    for _ in range(count):
        out.append(("random", perms_sample(), perms_sample()))
    # Grassmannian: odd-even interleave (shape staircase-ish), squared and against random
    g = list(range(1, n + 1, 2)) + list(range(2, n + 1, 2))
    out.append(("grass2", g, g))
    for _ in range(max(1, count // 4)):
        out.append(("grass_rand", g, perms_sample()))
    # many-descent, small code: w0 of S_n with blocks reversed, e.g. n ... pattern 6 2 1 5 4 3 generalized
    w = [n] + list(range(n // 2, 0, -1)) + list(range(n - 1, n // 2, -1))
    assert sorted(w) == list(range(1, n + 1)), w
    out.append(("desc2", w, w))
    for _ in range(max(1, count // 4)):
        out.append(("desc_rand", w, perms_sample()))
    # longest element: single monomial
    w0 = list(range(n, 0, -1))
    out.append(("w0_w0", w0, w0))
    return out


def gen():
    rng = random.Random(7)
    for n in range(6, 13):
        for fam, u, v in families(n, rng, count=40 if n <= 10 else 12):
            print(" ".join(map(str, u)) + " - " + " ".join(map(str, v)) + f"  # n={n} {fam}")


def summarize():
    rows = []
    for line in sys.stdin:
        parts = line.split()
        rs, cold, warm, nt, nt2 = map(int, parts[:5])
        tag = line.split("#")[-1].split()
        n = int(tag[0].split("=")[1]); fam = tag[1]
        rows.append((n, fam, rs, cold, warm, nt))
        if nt != nt2:
            print("MISMATCH", line.strip())
    fams = sorted(set(r[1] for r in rows), key=lambda f: ["random", "grass2", "grass_rand", "desc2", "desc_rand", "w0_w0"].index(f))
    print(f"{'family':>11} {'n':>3} {'count':>5} {'rs med us':>10} {'sm cold us':>11} {'sm warm us':>11} {'rs/cold':>8} {'rs/warm':>8} {'rs total s':>10} {'cold tot s':>10} {'nterms med':>10}")
    for fam in fams:
        for n in sorted(set(r[0] for r in rows)):
            sel = [r for r in rows if r[0] == n and r[1] == fam]
            if not sel:
                continue
            rs = np.array([r[2] for r in sel], float); cold = np.array([r[3] for r in sel], float); warm = np.array([r[4] for r in sel], float); nt = np.array([r[5] for r in sel])
            print(f"{fam:>11} {n:>3} {len(sel):>5} {np.median(rs)/1e3:10.1f} {np.median(cold)/1e3:11.1f} {np.median(warm)/1e3:11.1f} {np.exp(np.mean(np.log(rs/cold))):8.2f} {np.exp(np.mean(np.log(rs/warm))):8.2f} {rs.sum()/1e9:10.3f} {cold.sum()/1e9:10.3f} {int(np.median(nt)):>10}")


if __name__ == "__main__":
    gen() if sys.argv[1] == "gen" else summarize()
