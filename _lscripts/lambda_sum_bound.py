"""Bound Lambda(dom(u,v)) in terms of Lambda(dom u) and Lambda(dom v).

Lambda(w) = n! / prod_i d_i(w),  d_i(w) = 1 + #{j<i : w(j)<w(i)}.
For dominant sigma this equals |[e,sigma]_R|.

The definition in the paper is code(dom(u,v)) = code(dom u) + code(dom v),
so the whole question is about partitions: is
    Lambda(alpha + beta) <= Lambda(alpha) * Lambda(beta) ?
"""

import math
from itertools import permutations


def uncode(cd, N):
    cd = list(cd) + [0] * (N - len(cd))
    avail = list(range(1, N + 1))
    return [avail.pop(c) for c in cd]


def dvec(w):
    return [1 + sum(1 for j in range(i) if w[j] < w[i]) for i in range(len(w))]


def Lam_perm(w):
    p = 1
    for d in dvec(w):
        p *= d
    num = math.factorial(len(w))
    assert num % p == 0, (w, p)
    return num // p


def ambient(lam):
    return max([l + i + 1 for i, l in enumerate(lam)] + [1])


def Lam(lam):
    lam = [x for x in lam if x > 0]
    if not lam:
        return 1
    return Lam_perm(uncode(lam, ambient(lam)))


def inv_set(w):
    n = len(w)
    return frozenset((i, j) for i in range(n) for j in range(i + 1, n) if w[i] > w[j])


def right_interval_size(sigma):
    """|[e,sigma]_R| by brute force."""
    n = len(sigma)
    I = inv_set(sigma)
    return sum(1 for w in permutations(range(1, n + 1)) if inv_set(w) <= I)


def parts(maxpart, maxlen):
    """All partitions with at most maxlen parts, each at most maxpart."""
    out = []

    def rec(cur, prev, k):
        out.append(tuple(cur))
        if k == 0:
            return
        for v in range(1, min(prev, maxpart) + 1):
            rec(cur + [v], v, k - 1)

    rec([], maxpart, maxlen)
    return sorted(set(out))


def padding_stable():
    ok = True
    for lam in parts(3, 3):
        vals = set()
        for extra in range(4):
            N = ambient(lam) + extra
            vals.add(Lam_perm(uncode(lam, N)))
        if len(vals) != 1:
            ok = False
            print(f"  NOT STABLE {lam}: {vals}")
    print(f"padding stability of Lambda: {'OK' if ok else 'FAILED'}")


def matches_right_interval():
    bad = 0
    for lam in parts(3, 3):
        N = ambient(lam)
        if N > 7:
            continue
        s = uncode(lam, N)
        if Lam_perm(s) != right_interval_size(s):
            bad += 1
            print(f"  mismatch lam={lam} sigma={s} Lam={Lam_perm(s)} |[e,s]_R|={right_interval_size(s)}")
    print(f"Lambda == |[e,sigma]_R| on dominants: {bad} mismatches")


def multiset_formula():
    """Is Lambda(lam) = N!/prod_c m_c! for the multiset of parts padded to N?"""
    bad = []
    for lam in parts(3, 4):
        N = ambient(lam)
        padded = list(lam) + [0] * (N - len(lam))
        mult = {}
        for x in padded:
            mult[x] = mult.get(x, 0) + 1
        pred = math.factorial(N)
        for m in mult.values():
            pred //= math.factorial(m)
        if pred != Lam(lam):
            bad.append((lam, Lam(lam), pred))
    print(f"multiset/multinomial formula: {len(bad)} mismatches out of {len(parts(3, 4))}")
    for b in bad[:6]:
        print(f"   lam={b[0]}  Lambda={b[1]}  multinomial={b[2]}")


def product_conjecture(maxpart, maxlen):
    P = parts(maxpart, maxlen)
    worst_up = None
    worst_lo = None
    bad_up = bad_lo = 0
    for a in P:
        for b in P:
            k = max(len(a), len(b))
            aa = list(a) + [0] * (k - len(a))
            bb = list(b) + [0] * (k - len(b))
            s = tuple(x + y for x, y in zip(aa, bb))
            L, La, Lb = Lam(s), Lam(a), Lam(b)
            r = L / (La * Lb)
            if r > 1 + 1e-12:
                bad_up += 1
            if worst_up is None or r > worst_up[0]:
                worst_up = (r, a, b, L, La, Lb)
            r2 = max(La, Lb) / L
            if r2 > 1 + 1e-12:
                bad_lo += 1
            if worst_lo is None or r2 > worst_lo[0]:
                worst_lo = (r2, a, b, L, La, Lb)
    print(f"\npartitions <= {maxpart} with <= {maxlen} parts: {len(P)}, pairs {len(P) ** 2}")
    print(f"  Lambda(a+b) <= Lambda(a)Lambda(b):        {bad_up} failures")
    print(f"  max(Lambda(a),Lambda(b)) <= Lambda(a+b):  {bad_lo} failures")
    r, a, b, L, La, Lb = worst_up
    print(f"  tightest upper: a={a} b={b}  Lam(a+b)={L} Lam(a)Lam(b)={La * Lb}  ratio={r:.4f}")
    r, a, b, L, La, Lb = worst_lo
    print(f"  tightest lower: a={a} b={b}  Lam(a+b)={L} max={max(La, Lb)}  ratio={r:.4f}")


if __name__ == "__main__":
    padding_stable()
    matches_right_interval()
    multiset_formula()
    print("\nsmall values of Lambda:")
    for lam in parts(3, 3)[:20]:
        print(f"   lam={lam!s:16s} Lambda={Lam(lam)}")
    product_conjecture(3, 3)
    product_conjecture(4, 4)
