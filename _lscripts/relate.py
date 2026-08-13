"""Relationship between A_i,B_i,C_i in common S_M (M=2n-1).
A_i=sigma_a(i)-a_i, B_i=sigma_b(i)-b_i, C_i=sigma_g(i)-(a_i+b_i), g=a+b.
Identity: C_i = A_i + B_i - Delta_i,  Delta_i = sigma_a(i)+sigma_b(i)-sigma_g(i).

Test, exhaustively over partitions:
  (R1) Delta_i >= 0  ?   (i.e. C_i <= A_i + B_i pointwise)
  (R2) C_i >= 1 always (it must, it's a d-value)
  (R3) min(A_i,B_i) <= C_i <= max(A_i,B_i)  ?
  (R4) C_i >= A_i + B_i - i  ? (since sigma_g(i)<=... crude)
  (R5) sum Delta_i = M(M+1)/2 (tautology check)
"""
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


def dvec(code, M):
    ol = sigma_oneline(code, M)
    return [ol[i]-(code[i] if i < len(code) else 0) for i in range(M)], ol


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def main():
    maxpart, maxlen = 4, 4
    parts = [t for t in iproduct(range(maxpart+1), repeat=maxlen) if is_part(t)]
    n = maxlen + maxpart + 1
    M = 2*n - 1

    r1_fail = r3lo_fail = r3hi_fail = r4_fail = ident_fail = 0
    checks = 0
    r3lo_ex = r3hi_ex = r1_ex = None
    for a in parts:
        for b in parts:
            g = tuple(a[i]+b[i] for i in range(maxlen))
            A, sa = dvec(a, M); B, sb = dvec(b, M); C, sg = dvec(g, M)
            for i in range(M):
                checks += 1
                if C[i] != A[i]+B[i]-(sa[i]+sb[i]-sg[i]):
                    ident_fail += 1
                delta = sa[i]+sb[i]-sg[i]
                if delta < 0:
                    r1_fail += 1
                    if r1_ex is None: r1_ex = (a, b, i, A[i], B[i], C[i], delta)
                if C[i] < min(A[i], B[i]):
                    r3lo_fail += 1
                    if r3lo_ex is None: r3lo_ex = (a, b, i, A[i], B[i], C[i])
                if C[i] > max(A[i], B[i]):
                    r3hi_fail += 1
                    if r3hi_ex is None: r3hi_ex = (a, b, i, A[i], B[i], C[i])
                if C[i] < A[i]+B[i]-(i+1):
                    r4_fail += 1
    print(f"checks={checks}")
    print(f"(ident) C_i=A_i+B_i-Delta_i:        {ident_fail} fail")
    print(f"(R1) Delta_i>=0  (C_i<=A_i+B_i):     {r1_fail} fail" + (f"  ex a={r1_ex[0]} b={r1_ex[1]} i={r1_ex[2]} A={r1_ex[3]} B={r1_ex[4]} C={r1_ex[5]} d={r1_ex[6]}" if r1_ex else ""))
    print(f"(R3lo) C_i>=min(A_i,B_i):            {r3lo_fail} fail" + (f"  ex a={r3lo_ex[0]} b={r3lo_ex[1]} i={r3lo_ex[2]} A={r3lo_ex[3]} B={r3lo_ex[4]} C={r3lo_ex[5]}" if r3lo_ex else ""))
    print(f"(R3hi) C_i<=max(A_i,B_i):            {r3hi_fail} fail" + (f"  ex a={r3hi_ex[0]} b={r3hi_ex[1]} i={r3hi_ex[2]} A={r3hi_ex[3]} B={r3hi_ex[4]} C={r3hi_ex[5]}" if r3hi_ex else ""))
    print(f"(R4) C_i>=A_i+B_i-i:                 {r4_fail} fail")


if __name__ == "__main__":
    main()
