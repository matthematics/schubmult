"""The telescoping only ever compares bases differing by a PARTITION alpha.
Restricted claim (the one actually needed):

   R_i(alpha + S)  <=  R_i(S)      for all partitions alpha, S,
   whenever (alpha+S)+e_i and S+e_i are partitions,
   where R_i(gamma) = Lambda(sigma_{gamma+e_i})/Lambda(sigma_gamma).

Equivalently: adding a box in row i is cheaper on the heavier (by a partition) base.
This is exactly what makes Lambda(alpha+beta) <= Lambda(alpha)Lambda(beta).
"""
from math import factorial, comb
from itertools import product as iproduct
from fractions import Fraction
from schubmult import uncode


def Lam(code):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    L = len(code)
    if L == 0:
        return 1
    n = L + code[0] + 1
    perm = uncode(code)
    ol = list(perm) + list(range(len(list(perm)) + 1, n + 1))
    c = list(perm.code) + [0] * (n - len(perm.code))
    return factorial(n) // (lambda d: d)(_prod(ol[i] - c[i] for i in range(n)))


def _prod(it):
    p = 1
    for x in it:
        p *= x
    return p


def is_part(t):
    return all(t[i] >= t[i + 1] for i in range(len(t) - 1)) and all(x >= 0 for x in t)


def parts_upto(maxpart, maxlen):
    return [t for t in iproduct(range(maxpart + 1), repeat=maxlen) if is_part(t)]


def Ri(gamma, i):
    g2 = list(gamma); g2[i] += 1
    if not is_part(g2):
        return None
    return Fraction(Lam(g2), Lam(gamma))


def main():
    maxpart, maxlen = 5, 5
    parts = parts_upto(maxpart, maxlen)
    bad = 0; checks = 0; worst = None
    for alpha in parts:
        for S in parts:
            g = tuple(alpha[i] + S[i] for i in range(maxlen))
            if not is_part(g):  # alpha,S partitions => sum is a partition, always true
                continue
            for i in range(maxlen):
                rg = Ri(g, i)      # heavier base alpha+S
                rs = Ri(S, i)      # lighter base S
                if rg is None or rs is None:
                    continue
                checks += 1
                if rg > rs:        # want heavier base => smaller (or equal) marginal
                    bad += 1
                    if worst is None or (rg - rs) > worst[0]:
                        worst = (rg - rs, alpha, S, i, rg, rs)
    print(f"restricted antitone (base differs by partition alpha): {bad} failures over {checks} comparisons")
    if worst:
        d, alpha, S, i, rg, rs = worst
        print(f"  worst: i={i} alpha={alpha} S={S}  R_i(alpha+S)={float(rg):.4f} > R_i(S)={float(rs):.4f}")
    else:
        print("  HOLDS with no exceptions.")


if __name__ == "__main__":
    main()
