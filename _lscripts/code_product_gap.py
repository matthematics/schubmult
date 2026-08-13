"""How loose is  Lambda(sigma) <= prod_i (code_i(sigma)+1)  for dominant sigma?

Lambda(sigma) = |[e,sigma]_R| = # linear extensions of the forest P_sigma.
Compare with the Lehmer-code box product.  Equality iff ... ?
"""
import math, sys
from collections import defaultdict


def uncode(cd, N):
    cd = list(cd) + [0] * (N - len(cd))
    avail = list(range(1, N + 1))
    return [avail.pop(c) for c in cd]


def Lam_from_code(lam, N):
    s = uncode(lam, N)
    p = 1
    for i in range(N):
        p *= 1 + sum(1 for j in range(i) if s[j] < s[i])
    return math.factorial(N) // p


def parts_upto(N):
    """partitions lam with lam_i <= N-i (valid dominant codes in S_N)"""
    out = []
    def rec(cur, i, prev):
        out.append(tuple(cur))
        if i >= N:
            return
        for v in range(1, min(prev, N - 1 - i) + 1):
            rec(cur + [v], i + 1, v)
    rec([], 0, N)
    return out


def main(N):
    P = parts_upto(N)
    print(f"S_{N}: {len(P)} dominant permutations")
    worst = None
    eq_strict = eq_all = strict_total = 0
    byratio = defaultdict(int)
    for lam in P:
        L = Lam_from_code(lam, N)
        box = 1
        for x in lam:
            box *= x + 1
        assert L <= box, (lam, L, box)
        r = box / L
        if worst is None or r > worst[0]:
            worst = (r, lam, L, box)
        nz = [x for x in lam if x > 0]
        strict = all(nz[i] > nz[i + 1] for i in range(len(nz) - 1))
        if strict:
            strict_total += 1
            if L == box:
                eq_strict += 1
        if L == box:
            eq_all += 1
        byratio[round(r, 3)] += 1
    r, lam, L, box = worst
    print(f"  equality Lambda == prod(code+1): {eq_all} of {len(P)}")
    print(f"  strictly decreasing nonzero parts: {strict_total}, of which equality: {eq_strict}")
    print(f"  worst ratio box/Lambda = {r:.4f} at lam={lam}  (Lambda={L}, box={box})")
    tops = sorted(byratio.items(), key=lambda kv: -kv[0])[:6]
    print(f"  largest ratios seen: {tops}")
    # geometric mean of ratio
    import statistics
    rs = []
    for lam in P:
        L = Lam_from_code(lam, N)
        box = 1
        for x in lam:
            box *= x + 1
        rs.append(box / L)
    print(f"  mean ratio {statistics.mean(rs):.4f}, median {statistics.median(rs):.4f}")


if __name__ == "__main__":
    for N in range(2, (int(sys.argv[1]) if len(sys.argv) > 1 else 8) + 1):
        main(N)
