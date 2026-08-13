"""Everything padded into common S_M (M=2n-1). Print d-vectors A,B,C and look for
the arithmetic relation making prod A_i prod B_i <= M! prod C_i.

d_i(sigma_lambda) = sigma_lambda(i) - lambda_i = 1 + #{j<i : sigma(j)<sigma(i)}.
Padding with fixed points: a fixed point at position p contributes d_p = p.
"""
from math import factorial
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
    return [ol[i] - (code[i] if i < len(code) else 0) for i in range(M)]


def is_part(t):
    return all(t[i] >= t[i+1] for i in range(len(t)-1)) and all(x >= 0 for x in t)


def main():
    samples = [
        ((1,1),(1,1)),
        ((2,1),(1,1)),
        ((2,2),(1,1)),
        ((2,1),(2,1)),
        ((3,1),(1,1,1)),
        ((2,1,1),(1,1)),
        ((3,2,1),(2,1)),
    ]
    for a,b in samples:
        n = max(len(a),len(b)) + max(max(a),max(b)) + 1
        M = 2*n - 1
        A = dvec(a, M); B = dvec(b, M)
        g = tuple(a[i]+b[i] for i in range(max(len(a),len(b))))
        C = dvec(g, M)
        pA=1;pB=1;pC=1
        for i in range(M):
            pA*=A[i];pB*=B[i];pC*=C[i]
        rho = Fraction(pA*pB, factorial(M)*pC)
        print(f"alpha={a} beta={b} gamma={g}  M={M}")
        print(f"  A={A}")
        print(f"  B={B}")
        print(f"  C={C}")
        print(f"  A_i+B_i = {[A[i]+B[i] for i in range(M)]}")
        print(f"  A_i*B_i = {[A[i]*B[i] for i in range(M)]}")
        print(f"  C_i     = {C}")
        print(f"  sorted(A_i*B_i) = {sorted(A[i]*B[i] for i in range(M))}")
        print(f"  sorted(C_i)     = {sorted(C)}")
        print(f"  prodA*prodB={pA*pB}  M!*prodC={factorial(M)*pC}  rho={float(rho):.4f}")
        # Does C majorize / relate to A,B? check prod C_i vs prod max(A_i,B_i), prod (A_i+B_i-1)
        pmax=1;pAB1=1
        for i in range(M):
            pmax*=max(A[i],B[i]); pAB1*=(A[i]+B[i]-1)
        print(f"  prod max(A,B)={pmax}  prod(A+B-1)={pAB1}")
        print()


if __name__=="__main__":
    main()
