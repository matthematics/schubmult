"""Focused test of the denominator claim under the restriction ell(v) >= n.

claim2:  for every simple reflection s_k and every v with ell(v) >= n,
         N(s_k, v) >= n - 1.

We also bucket the MINIMAL N(s_k,v) by ell(v) to locate the exact threshold,
and report the strongest per-step inductive quantity actually needed.

N(u,v) = | U_{v'<=_L v} supp(S_u S_{v'}) |.
"""
import sys
from schubmult import Permutation
from schubmult.rings.schubert.schubert_ring import Sx


def left_downset(perm, n):
    start = tuple(perm[i] for i in range(n))
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        pos = [0] * (n + 2)
        for idx, val in enumerate(a):
            pos[val] = idx
        for i in range(1, n):
            if pos[i] > pos[i + 1]:
                b = list(a)
                b[pos[i]] = i + 1
                b[pos[i + 1]] = i
                b = tuple(b)
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
    return seen


_supp = {}


def supp(u, w):
    key = (u, w)
    r = _supp.get(key)
    if r is None:
        r = frozenset((Sx(Permutation(list(u))) * Sx(Permutation(list(w)))).keys())
        _supp[key] = r
    return r


def length(a, n):
    return sum(1 for i in range(n) for j in range(i + 1, n) if a[i] > a[j])


def run(n):
    perms = [tuple(p[i] for i in range(n)) for p in Permutation.all_permutations(n)]
    simple = []
    for k in range(n - 1):
        a = list(range(1, n + 1))
        a[k], a[k + 1] = a[k + 1], a[k]
        simple.append(tuple(a))

    # min N(s_k,v) bucketed by ell(v)
    by_len = {}   # ell -> (minN, example)
    for v in perms:
        lv = length(v, n)
        ds = left_downset(Permutation(list(v)), n)
        for s in simple:
            acc = set()
            for w in ds:
                acc |= supp(s, tuple(w))
            nv = len(acc)
            cur = by_len.get(lv)
            if cur is None or nv < cur[0]:
                by_len[lv] = (nv, (s, v))

    print(f"===== n={n} =====")
    print(f"  ell(v):  min_over_(s_k) N(s_k,v)   (need >= n-1 = {n-1})")
    threshold = None
    for lv in sorted(by_len):
        mn, ex = by_len[lv]
        ok = mn >= n - 1
        if ok and threshold is None:
            # first ell where it holds; but we want the LAST ell where it FAILS
            pass
        flag = "OK " if ok else "FAIL"
        print(f"    ell={lv:2d}: min N = {mn:3d}  [{flag}]  e.g. s={ex[0]} v={ex[1]}")
    # find max ell that still fails
    fails = [lv for lv, (mn, _) in by_len.items() if mn < n - 1]
    if fails:
        print(f"  ==> claim2 FAILS for ell(v) in {sorted(fails)}; holds for ell(v) >= {max(fails)+1}")
    else:
        print(f"  ==> claim2 holds for ALL v")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["5", "6"])):
        run(n)
