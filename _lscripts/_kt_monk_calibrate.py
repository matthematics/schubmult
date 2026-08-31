"""Scratch: calibrate Lenart-Postnikov Cor 8.2 against DoubleGrothendieckRing."""

import itertools

from sympy import cancel, expand, simplify

from schubmult import Permutation, uncode
from schubmult.rings.schubert.double_grothendieck_ring import DoubleGrothendieckRing
from schubmult.symbolic import sympify_sympy
from schubmult.symbolic.poly.variables import GeneratingSet

x = GeneratingSet("x")
y = GeneratingSet("y")
R = DoubleGrothendieckRing(x, y)
beta = R.beta


def chain(i, n):
    """(-omega_i)-chain of reflections: (a,b) with a<=i<b, a increasing, b decreasing."""
    return [(a, b) for a in range(1, i + 1) for b in range(n, i, -1)]


def monk_terms(v, i, n):
    """sum_J (-1)^{|J|-1} [w(J)] as a dict {perm: int}, J over ALL subsets."""
    L = chain(i, n)
    out = {}

    def rec(pos, w, size):
        if pos == len(L):
            out[w] = out.get(w, 0) + (-1) ** (size - 1)
            return
        rec(pos + 1, w, size)
        a, b = L[pos]
        w2 = w.swap(a - 1, b - 1)
        if w2.inv == w.inv + 1:
            rec(pos + 1, w2, size + 1)

    rec(0, Permutation(v), 0)
    return {k: c for k, c in out.items() if c}


def truth(v, i):
    return R(uncode([0] * (i - 1) + [1])) * R(v)


def report(v, i, n):
    v = Permutation(v)
    t = truth(v, i)
    terms = monk_terms(v, i, n)
    print(f"v={list(v)} i={i} n={n}")
    print("  truth   :", {tuple(k): simplify(sympify_sympy(c)) for k, c in t.items()})
    print("  combin  :", {tuple(k): c for k, c in terms.items()})
    # solve for C from a term w != v
    for w, c in terms.items():
        if w != v:
            got = t.get(w, 0)
            print(f"  C from {tuple(w)}: {cancel(sympify_sympy(got) / c)}")
            break


if __name__ == "__main__":
    for n in (3, 4):
        for v in Permutation.all_permutations(n):
            for i in range(1, n):
                report(v, i, n)
