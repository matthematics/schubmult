from schubmult import *
from schubmult.mult.groth_double import *
from schubmult.symbolic import *
import itertools

from sympy import Poly, fraction, together


def sanify(v, subs_dict):
    """To s/t space: (num, den) with num a genuine s/t polynomial, den free of s/t."""
    e = together(sympify_sympy(efficient_subs(v, subs_dict)))
    num, den = fraction(e)
    if any(str(g).startswith(("s_", "t_")) for g in den.free_symbols):
        raise ValueError("denominator not cleared")
    return num.expand(), den


def desanify(num, den, y, z, beta):
    """Monomial-wise back-map s_i -> y_i/(1+b y_i), t_j -> 1+b z_j; no cancel/GCD."""
    st = [g for g in num.free_symbols if str(g).startswith(("s_", "t_"))]
    if not st:
        return sympify(num) / sympify(den)
    p = Poly(num, *st)
    degs = {g: p.degree(g) for g in st}
    extra_den = S.One
    for g in st:
        nm = str(g)
        if nm.startswith("s_"):
            extra_den *= (S.One + beta * y[int(nm.split("_")[1])]) ** degs[g]
    total = S.Zero
    for monom, coeff in p.terms():
        term = sympify(coeff)
        for g, a in zip(p.gens, monom):
            nm = str(g)
            idx = int(nm.split("_")[1])
            if nm.startswith("s_"):
                term *= y[idx] ** a * (S.One + beta * y[idx]) ** (degs[g] - a)
            else:
                term *= (S.One + beta * z[idx]) ** a
        total += term
    return total / (sympify(den) * extra_den)


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    t = GeneratingSet("t")
    y = GeneratingSet("y")
    z = GeneratingSet("z")
    s = GeneratingSet("s")
    beta = Gx._beta
    subs_dict = {z[i]: (t[i] - S.One)/beta for i in range(50)}
    subs_dict |= {y[i]: s[i]/(S.One - beta * s[i]) for i in range(50)}
    for perm1, perm2 in itertools.product(perms, repeat=2):
        prd = DGx(perm1) * DGx(perm2, "z")
        prd2 = DGx([]).ring.from_dict({k: desanify(*sanify(v, subs_dict), y, z, beta) for k, v in prd.items()})
        print(prd2)