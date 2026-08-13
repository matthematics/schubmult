"""Lambda(alpha+beta) <= Lambda(alpha)Lambda(beta) for partitions (code sum).

Also tests a candidate proof route:
  (R1) h_i(a+b) == min(h_i(a), h_i(b))  on a common padding
  (R2) prod_i max(h_i(a),h_i(b)) <= N!
(R1)+(R2) together imply the inequality, since
  Lam(a)Lam(b)/Lam(a+b) = prod_i h_i(a+b)*N!/(h_i(a)h_i(b)) ... rearranged.
"""

import math
from itertools import permutations


def uncode(cd, N):
    cd = list(cd) + [0] * (N - len(cd))
    avail = list(range(1, N + 1))
    return [avail.pop(c) for c in cd]


def hvec(lam, N):
    s = uncode(lam, N)
    return [1 + sum(1 for j in range(i) if s[j] < s[i]) for i in range(N)]


def Lam(lam, N):
    p = 1
    for h in hvec(lam, N):
        p *= h
    return math.factorial(N) // p


def ambient(lam):
    return max([l + i + 1 for i, l in enumerate(lam)] + [1])


def parts(maxpart, maxlen):
    out = []

    def rec(cur, prev, k):
        out.append(tuple(cur))
        if k:
            for v in range(1, min(prev, maxpart) + 1):
                rec(cur + [v], v, k - 1)

    rec([], maxpart, maxlen)
    return sorted(set(out))


def run(maxpart, maxlen, verbose=False):
    P = parts(maxpart, maxlen)
    bad = badR1 = badR2 = 0
    worst = None
    tot = 0
    for a in P:
        for b in P:
            k = max(len(a), len(b))
            aa = list(a) + [0] * (k - len(a))
            bb = list(b) + [0] * (k - len(b))
            s = tuple(x + y for x, y in zip(aa, bb))
            N = max(ambient(s), ambient(a), ambient(b))
            tot += 1
            La, Lb, Ls = Lam(a, N), Lam(b, N), Lam(s, N)
            if Ls > La * Lb:
                bad += 1
                if verbose and bad <= 3:
                    print(f"   FAIL a={a} b={b}: {Ls} > {La}*{Lb}")
            if a and b:
                r = Ls / (La * Lb)
                if worst is None or r > worst[0]:
                    worst = (r, a, b, Ls, La * Lb)
            ha, hb, hs = hvec(aa, N), hvec(bb, N), hvec(s, N)
            if hs != [min(x, y) for x, y in zip(ha, hb)]:
                badR1 += 1
                if verbose and badR1 <= 3:
                    print(f"   R1 fail a={a} b={b} N={N}: h(a)={ha} h(b)={hb} h(a+b)={hs}")
            pm = 1
            for x, y in zip(ha, hb):
                pm *= max(x, y)
            if pm > math.factorial(N):
                badR2 += 1
    print(f"parts<={maxpart}, len<={maxlen}: {tot} pairs")
    print(f"   Lambda(a+b) <= Lambda(a)Lambda(b) : {bad} failures")
    print(f"   R1  h(a+b)==min(h(a),h(b))        : {badR1} failures")
    print(f"   R2  prod max(h_i) <= N!           : {badR2} failures")
    if worst:
        r, a, b, Ls, pr = worst
        print(f"   tightest nontrivial: a={a} b={b}  {Ls} vs {pr}  ratio={r:.4f}")


if __name__ == "__main__":
    run(3, 3, verbose=True)
    run(4, 4)
    run(5, 5)
    run(6, 2)
    run(2, 6)
