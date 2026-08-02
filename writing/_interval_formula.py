import math
from functools import lru_cache
from itertools import permutations

from schubmult import Permutation, uncode


def dominant_from_code(nu, n):
    return uncode(list(nu) + [0] * (n - len(nu)))


def interval_size(mu, n):
    """|[e,mu]_L| = # linear extensions of the non-inversion poset of mu^{-1}."""
    sigma = ~mu
    s = [sigma[i] for i in range(n)]
    below = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if s[i] < s[j]:
                below[j] |= 1 << i
    for j in range(n):
        changed = True
        while changed:
            changed = False
            for i in range(n):
                if (below[j] >> i) & 1 and (below[i] & ~below[j]):
                    below[j] |= below[i]
                    changed = True
    below = tuple(below)
    full = (1 << n) - 1

    @lru_cache(maxsize=None)
    def f(mask):
        if mask == full:
            return 1
        tot = 0
        for j in range(n):
            if not (mask >> j) & 1 and (below[j] & ~mask) == 0:
                tot += f(mask | (1 << j))
        return tot

    r = f(0)
    f.cache_clear()
    return r


def blocks(nu):
    out = []
    for x in [y for y in nu if y != 0]:
        if out and out[-1][0] == x:
            out[-1][1] += 1
        else:
            out.append([x, 1])
    return out


def closed_form(nu):
    return math.prod(math.comb(k + m, m) for k, m in blocks(nu))


def prodplus1(nu):
    return math.prod(x + 1 for x in nu if x != 0)


def conjugate(nu):
    nu = [x for x in nu if x != 0]
    return [sum(1 for x in nu if x > j) for j in range(nu[0])] if nu else []


def partitions_upto(maxpart, maxlen):
    def rec(mx, ln):
        yield []
        if ln:
            for p in range(1, mx + 1):
                for rest in rec(p, ln - 1):
                    yield [p, *rest]

    seen = set()
    for p in rec(maxpart, maxlen):
        if tuple(p) not in seen:
            seen.add(tuple(p))
            yield p


print("=== 1. linear-extension count == direct enumeration of [e,mu]_L ===")
bad = 0
for n in range(1, 8):
    allp = [Permutation(list(p)) for p in permutations(range(1, n + 1))]
    for mu in allp:
        if mu.is_dominant and sum(1 for w in allp if (w * (~mu)).inv == mu.inv - w.inv) != interval_size(mu, n):
            bad += 1
print("  mismatches:", bad)

print("\n=== 2. closed form  prod_blocks C(k+m,m)  on code(mu^{-1}) ===")
badc = badu = tight = tot = 0
for nu in partitions_upto(6, 6):
    n = sum(nu) + len(nu) + 2
    if not nu or n > 13:
        continue
    mu = dominant_from_code(nu, n)
    size = interval_size(mu, n)
    nusig = list((~mu).trimcode)
    tot += 1
    if closed_form(nusig) != size:
        badc += 1
        if badc <= 5:
            print("   MISS", nu, "-> code(mu^-1)", nusig, "true", size, "pred", closed_form(nusig))
    if prodplus1(nusig) < size:
        badu += 1
    if prodplus1(nusig) == size:
        tight += 1
print(f"  partitions tested: {tot}   closed-form mismatches: {badc}")
print(f"  prod(nu_i+1) fails as upper bound: {badu}    exactly tight in {tight}")

print("\n=== 3. code(mu^{-1}) == conjugate(code(mu)) for dominant mu? ===")
bad = sum(
    1
    for nu in partitions_upto(6, 6)
    if nu and sum(nu) + len(nu) + 2 <= 13 and list((~dominant_from_code(nu, sum(nu) + len(nu) + 2)).trimcode) != conjugate(nu)
)
print("  mismatches:", bad)

print("\n=== 4. closed form invariant under conjugation? ===")
bad = sum(1 for nu in partitions_upto(7, 7) if nu and closed_form(nu) != closed_form(conjugate(nu)))
print("  mismatches:", bad)

print("\n=== 5. doubled staircase (2n-2, 2n-4, ..., 2) ===")
for n in range(2, 9):
    nu = [2 * (n - i) - 2 for i in range(1, n)]
    print(f"  n={n} code={nu} closed={closed_form(nu)} (2n-1)!!={math.prod(range(1, 2 * n, 2))} prod+1={prodplus1(nu)}")
