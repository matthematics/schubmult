"""Verify the clean proof steps for the Lemma  Lambda(mu_v) <= n^L N(u,v).

Key claims to verify exhaustively (n=4,5,6):
  (A)  w <=_R delta  =>  code(w)_i <= code(delta)_i  for all i   (delta any perm)
       => |[e,delta]_R| <= prod_i(code(delta)_i + 1).
  (C)  Lambda(mu_v) = |[e,mu_v]_L| = |[e,dom(v^{-1})]_R|.
  (D)  |[e,dom(v^{-1})]_R| <= prod_i(theta_i+1) <= n^L,   theta=theta(v^{-1}),
       L = #nonzero parts of theta.
  (E)  the descent-count version prod over DISTINCT parts (n^k) can be exceeded,
       i.e. Lambda(mu_v) > n^k for some v  (shows k is wrong, L is right).
"""
import sys
from math import prod
from schubmult import Permutation


def right_downset(perm, n):
    start = tuple(perm[i] for i in range(n))
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        for i in range(n - 1):
            if a[i] < a[i + 1]:  # right ascent -> can add inversion (go up); we go DOWN: swap a descent
                pass
        # generate right-weak-order DOWN-covers: swap adjacent positions that ARE a descent
        for i in range(n - 1):
            if a[i] > a[i + 1]:
                b = list(a)
                b[i], b[i + 1] = b[i + 1], b[i]
                b = tuple(b)
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
    return seen


def code_of(a, n):
    return [sum(1 for j in range(i + 1, n) if a[j] < a[i]) for i in range(n)]


def run(n):
    okA = okD = True
    max_ratio_nL = 0.0
    fail_k = None
    for v in Permutation.all_permutations(n):
        vi = ~v
        theta = sorted((c for c in vi.code), reverse=True)
        L = len([c for c in theta if c != 0])
        k = len(set(c for c in theta if c != 0))  # distinct nonzero parts = #descents of dom
        delta = vi.minimal_dominant_above()  # = dom(v^{-1}), dominant
        da = tuple(delta[i] for i in range(n))
        cdelta = code_of(da, n)
        ds = right_downset(delta, n)
        # (A) code componentwise for every w <=_R delta
        for w in ds:
            cw = code_of(w, n)
            if any(cw[i] > cdelta[i] for i in range(n)):
                okA = False
        Lambda = len(ds)
        prodbound = prod(c + 1 for c in cdelta)
        nL = n ** L
        nk = n ** k
        # (D)
        if not (Lambda <= prodbound <= nL):
            okD = False
        max_ratio_nL = max(max_ratio_nL, Lambda / nL)
        # (E) does n^k fail?
        if Lambda > nk and fail_k is None:
            fail_k = (tuple(v), tuple(da), theta, L, k, Lambda, nk)
    print(f"n={n}: (A) code-componentwise = {okA};  (D) Lambda<=prod(c+1)<=n^L = {okD}")
    print(f"      max Lambda/n^L = {max_ratio_nL:.4f}   (<=1 confirms Lambda<=n^L)")
    print(f"      n^k (distinct-parts) failure example: {fail_k}")


if __name__ == "__main__":
    for n in (int(a) for a in (sys.argv[1:] or ["4", "5", "6"])):
        run(n)
