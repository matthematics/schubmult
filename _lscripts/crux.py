"""THE CRUX for an elementary proof of rho<=1:

Pad to common S_M. Claim:  sigma_alpha <=_R sigma_{alpha+beta}  AND  sigma_beta <=_R sigma_{alpha+beta}
i.e. Inv(sigma_alpha) subset Inv(sigma_gamma), Inv(sigma_beta) subset Inv(sigma_gamma),
for all partitions alpha,beta with gamma=alpha+beta.

If TRUE then E_gamma superset E_alpha cup E_beta, and the union bound
   P[E_gamma] >= P[E_alpha cup E_beta] >= P[E_alpha]P[E_beta]
gives Lambda(alpha)Lambda(beta) <= M! Lambda(gamma) = rho<=1.  COMPLETE PROOF.

Also directly re-confirm the final rho<=1 to be safe.
"""
from math import factorial
from itertools import product as iproduct
from schubmult import uncode


def sigma_oneline(code, M):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        return list(range(1, M+1))
    ol = list(uncode(code))
    return ol + list(range(len(ol)+1, M+1))


def inv_set(ol):
    M = len(ol)
    return frozenset((i, j) for i in range(M) for j in range(i+1, M) if ol[i] > ol[j])


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def Lam(code, M):
    ol = sigma_oneline(code, M)
    p = 1
    for i in range(M):
        p *= ol[i] - (code[i] if i < len(code) else 0)
    return factorial(M)//p


def main():
    maxpart, maxlen = 5, 4
    parts = [t for t in iproduct(range(maxpart+1), repeat=maxlen) if is_part(t)]
    n = maxlen + maxpart + 1
    M = 2*n - 1
    fact = factorial(M)

    crux_a_fail = crux_b_fail = 0
    rho_fail = 0
    checks = 0
    ex = None
    for a in parts:
        Ia = inv_set(sigma_oneline(a, M))
        for b in parts:
            g = tuple(a[i]+b[i] for i in range(maxlen))
            Ib = inv_set(sigma_oneline(b, M))
            Ig = inv_set(sigma_oneline(g, M))
            checks += 1
            if not (Ia <= Ig):
                crux_a_fail += 1
                if ex is None: ex = ('a', a, b, g)
            if not (Ib <= Ig):
                crux_b_fail += 1
                if ex is None: ex = ('b', a, b, g)
            if Lam(a, M)*Lam(b, M) > fact*Lam(g, M):
                rho_fail += 1
    print(f"M={M}, {checks} pairs")
    print(f"  CRUX  sigma_alpha <=_R sigma_gamma:  {crux_a_fail} fail")
    print(f"  CRUX  sigma_beta  <=_R sigma_gamma:  {crux_b_fail} fail")
    print(f"  rho<=1 (Lam(a)Lam(b)<=M! Lam(g)):    {rho_fail} fail")
    if ex:
        print(f"  first crux failure: {ex}")
    else:
        print("  CRUX HOLDS with no exceptions -> elementary proof of rho<=1 is valid.")


if __name__ == "__main__":
    main()
