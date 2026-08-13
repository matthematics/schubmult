"""User's route: C_i >= sqrt(A_i B_i) pointwise would give rho<=1, since
   prod C_i >= sqrt(prod A_i prod B_i) >= prod A_i prod B_i / M!   (each prod<=M!).
Condition offered: Delta_i < (A_i+B_i)/2  => C_i=A+B-Delta > (A+B)/2 >= sqrt(AB).

TESTS:
 (P1) pointwise  C_i >= sqrt(A_i B_i)  ?           (equiv Delta_i <= (A_i+B_i)-sqrt(A_iB_i))
 (P1b) pointwise Delta_i <= (A_i+B_i)/2 ?          (the offered sufficient condition)
 (P2) PRODUCT   prod C_i >= sqrt(prod A_i prod B_i) ?   (all that's really needed)
 (P3) sanity    prod A_i <= M! and prod B_i <= M! ?
"""
from math import factorial, isqrt
from itertools import product as iproduct
from fractions import Fraction
from schubmult import uncode


def sigma_oneline(code, M):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        return list(range(1, M+1))
    ol = list(uncode(code))
    return ol + list(range(len(ol)+1, M+1))


def dvec(code, M):
    ol = sigma_oneline(code, M)
    return [ol[i]-(code[i] if i < len(code) else 0) for i in range(M)]


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def main():
    maxpart, maxlen = 4, 4
    parts = [t for t in iproduct(range(maxpart+1), repeat=maxlen) if is_part(t)]
    n = maxlen + maxpart + 1
    M = 2*n - 1
    fact = factorial(M)

    p1_fail = p1b_fail = p2_fail = p3_fail = 0
    checks = 0; prods = 0
    p1_ex = p2_ex = None
    for a in parts:
        for b in parts:
            g = tuple(a[i]+b[i] for i in range(maxlen))
            A = dvec(a, M); B = dvec(b, M); C = dvec(g, M)
            pA = pB = pC = 1
            for i in range(M):
                checks += 1
                # P1 pointwise sqrt: C_i^2 >= A_i B_i
                if C[i]*C[i] < A[i]*B[i]:
                    p1_fail += 1
                    if p1_ex is None: p1_ex = (a, b, i, A[i], B[i], C[i])
                # P1b Delta <= (A+B)/2  i.e. 2*(A+B-C) <= A+B  i.e. A+B <= 2C
                if A[i]+B[i] > 2*C[i]:
                    p1b_fail += 1
                pA *= A[i]; pB *= B[i]; pC *= C[i]
            prods += 1
            # P2 product sqrt: (prod C)^2 >= prod A prod B
            if pC*pC < pA*pB:
                p2_fail += 1
                if p2_ex is None: p2_ex = (a, b, Fraction(pC*pC, pA*pB))
            if pA > fact: p3_fail += 1
            if pB > fact: p3_fail += 1
    print(f"M={M}, {prods} pairs, {checks} pointwise checks")
    print(f"(P1) pointwise C_i>=sqrt(A_iB_i):      {p1_fail} fail" + (f"  ex a={p1_ex[0]} b={p1_ex[1]} i={p1_ex[2]} A={p1_ex[3]} B={p1_ex[4]} C={p1_ex[5]}" if p1_ex else ""))
    print(f"(P1b) pointwise Delta_i<=(A_i+B_i)/2:  {p1b_fail} fail")
    print(f"(P2) PRODUCT prodC>=sqrt(prodA prodB): {p2_fail} fail" + (f"  worst (prodC^2/prodAB)={float(p2_ex[2]):.4f} a={p2_ex[0]} b={p2_ex[1]}" if p2_ex else ""))
    print(f"(P3) prodA,prodB <= M!:                {p3_fail} fail")


if __name__ == "__main__":
    main()
