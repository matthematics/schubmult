"""Search for the correct exponent in the INTERVAL bound
   |[e,mu_v]_L| <= n^param * |[e,v]_L|
for several candidate parameters, to identify the right one.

Candidates:
  k     = #descents(v)              (right descents)
  kinv  = #descents(v^{-1})
  Lpos  = # nonzero parts of theta(v^{-1})   (= #nonzero code entries of v^{-1})
  Ldist = # distinct nonzero parts of theta(v^{-1}) (= #descents of dom(v^{-1}) = #descents(mu_v))
"""
import math
from schubmult import Permutation


def left_weak_downset_size(perm, n):
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
    return len(seen)


def run(n):
    best = {"k": (0.0, None), "kinv": (0.0, None), "Lpos": (0.0, None), "Ldist": (0.0, None), "D": (0.0, None)}
    maxestar = (0.0, None)
    d_lt_L = []  # cases where maxdesc(v) < Lpos  (would break n^D as an upper bound)
    for v in Permutation.all_permutations(n):
        vi = ~v
        k = len(v.descents())
        kinv = len(vi.descents())
        theta = [c for c in vi.theta()]
        Lpos = len([c for c in theta if c != 0])
        Ldist = len(set(c for c in theta if c != 0))
        desc = v.descents()
        D = (max(desc) + 1) if desc else 0  # max descent of v = #x-variables of S_v
        if D < Lpos:
            d_lt_L.append((tuple(v), D, Lpos))
        mu = ~(vi.minimal_dominant_above())
        lv = left_weak_downset_size(v, n)
        lmu = left_weak_downset_size(mu, n)
        estar = math.log(lmu / lv) / math.log(n) if lmu > lv else 0.0
        info = (tuple(v), tuple(mu), k, kinv, Lpos, Ldist, D, lv, lmu)
        if estar > maxestar[0]:
            maxestar = (estar, info)
        for name, p in (("k", k), ("kinv", kinv), ("Lpos", Lpos), ("Ldist", Ldist), ("D", D)):
            r = lmu / (n ** p * lv) if p > 0 else lmu / lv
            if r > best[name][0]:
                best[name] = (r, info)
    print(f"===== n={n} =====")
    print(f"  max needed exponent e* = {maxestar[0]:.4f}   {maxestar[1]}")
    print(f"  cases maxdesc(v) < Lpos: {len(d_lt_L)}"
          + (f"  e.g. {d_lt_L[:3]}" if d_lt_L else "  (maxdesc(v) >= Lpos always)"))
    for name in ("k", "kinv", "Lpos", "Ldist", "D"):
        r, info = best[name]
        flag = "OK " if r <= 1.0 + 1e-9 else "FAIL"
        print(f"  [{flag}] max |[e,mu]_L|/(n^{name} |[e,v]_L|) = {r:.4f}   {info}")


if __name__ == "__main__":
    import sys
    ns = [int(a) for a in sys.argv[1:]] or list(range(3, 8))
    for n in ns:
        run(n)
