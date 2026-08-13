"""Explore the forest (tree) structure of the dominant permutation of a partition code,
and how F_{alpha+beta} is built from F_alpha and F_beta.

Forest P_sigma (sigma dominant, 132-avoiding): on positions 1..N,
 i is BELOW j (j covers toward root) iff ... down-set of i = {j<=i : sigma(j)<=sigma(i)} is a chain.
parent(i) = the largest j < i with sigma(j) < sigma(i)  (immediate ancestor in the chain).
Root(s): nodes with no such j.
down-set size h_i = |{j<=i : sigma(j)<=sigma(i)}| = sigma(i) - lambda_i.
Lambda = N!/prod h_i.
"""
from math import factorial
from schubmult import uncode


def dominant_oneline(code):
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    if not code:
        return [1]
    perm = uncode(code)
    return list(perm)


def forest(code):
    ol = dominant_oneline(code)
    N = len(ol)
    sigma = ol
    parent = [None] * N
    h = [0] * N
    for i in range(N):
        below = [j for j in range(i) if sigma[j] < sigma[i]]  # j<i, sigma(j)<sigma(i)
        h[i] = 1 + len(below)
        parent[i] = (max(below) if below else None)
    return sigma, parent, h


def Lam(code):
    sigma, parent, h = forest(code)
    N = len(sigma)
    p = 1
    for x in h:
        p *= x
    return factorial(N) // p


def show(code, label):
    sigma, parent, h = forest(code)
    N = len(sigma)
    print(f"  {label}: code={tuple(code)} sigma(one-line)={sigma}")
    print(f"       parents={[ (p+1 if p is not None else 0) for p in parent]}  (0=root, 1-indexed)")
    print(f"       hooks  ={h}   Lambda={Lam(code)}  N={N}")


def main():
    samples = [
        ((1,1), (1,1)),
        ((2,1), (1,1)),
        ((2,2), (1,1)),
        ((2,1), (2,1)),
        ((3,1), (1,1,1)),
        ((2,1,1), (1,1)),
    ]
    for a, b in samples:
        L = max(len(a), len(b))
        aa = list(a)+[0]*(L-len(a)); bb=list(b)+[0]*(L-len(b))
        g = tuple(aa[i]+bb[i] for i in range(L))
        print(f"alpha={a} beta={b} -> gamma={g}")
        show(a, "F_alpha")
        show(b, "F_beta ")
        show(g, "F_gamma")
        # node counts
        Na = len(dominant_oneline(a)); Nb=len(dominant_oneline(b)); Ng=len(dominant_oneline(g))
        print(f"   node counts: N_a={Na} N_b={Nb} N_gamma={Ng}   (2n-1 style: N_a+N_b-? )")
        print(f"   rho = {Lam(g)}/({Lam(a)}*{Lam(b)}) = {Lam(g)/(Lam(a)*Lam(b)):.4f}")
        print()


if __name__ == "__main__":
    main()
