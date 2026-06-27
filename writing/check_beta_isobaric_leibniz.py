#!/usr/bin/env python3
"""Check the beta-isobaric Leibniz identity on Schubert polynomials.

Identity tested (for i >= 1):
    pi_i(f g) = pi_i(f) s_i(g) + f pi_i(g)
where
    pi_i = partial_i (1 + beta x_{i+1}).

This script uses the NilHecke beta-isobaric factory:
    NilHeckeRing.isobaric(..., groth=True, groth_beta=beta)
"""

from __future__ import annotations

from schubmult import Permutation, Sx
from schubmult.abc import x, beta
from schubmult.rings.schubert.nil_hecke import NilHeckeRing
from schubmult.symbolic import S, Symbol, expand


def simple_reflection(i: int) -> Permutation:
    """Return s_i in 1-based indexing."""
    if i < 1:
        raise ValueError("i must be >= 1")
    return Permutation([]).swap(i - 1, i)


def s_i_action(elem, i: int):
    """Apply s_i to a Schubert element by swapping x_i and x_{i+1}."""
    poly = elem.as_polynomial()
    swapped = poly.subs({x[i]: x[i + 1], x[i + 1]: x[i]})
    return elem.ring.from_expr(swapped)


nh = NilHeckeRing(x)
def beta_isobaric_operator(i: int, beta):
    """Build pi_i = partial_i(1 + beta x_{i+1}) as a NilHecke element."""
    si = simple_reflection(i)
    return nh.isobaric(si, groth=True, groth_beta=beta)


def schubert_basis_samples(max_n: int = 5):
    """Generate a small set of Schubert basis elements for testing."""
    out = []
    seen = set()
    for n in range(1, max_n + 1):
        for p in Permutation.all(n):
            if p in seen:
                continue
            seen.add(p)
            out.append(Sx(p))
    return out


def check_identity(i: int, beta, verbose: bool = False) -> tuple[int, int]:
    pi_i = beta_isobaric_operator(i, beta)
    elems = schubert_basis_samples(max_n=5)

    tested = 0
    passed = 0
    for f in elems:
        for g in elems:
            lhs = pi_i.apply(f * g)
            rhs = pi_i.apply(f) * s_i_action(g, i) + f * pi_i.apply(g)
            diff = expand(lhs.as_polynomial() - rhs.as_polynomial(), deep=True)
            tested += 1
            if diff == S.Zero:
                passed += 1
            elif verbose:
                print("FAILED")
                print("i=", i, "f=", f, "g=", g)
                print("lhs=", expand(lhs.as_polynomial(), deep=True))
                print("rhs=", expand(rhs.as_polynomial(), deep=True))
                print("diff=", diff)
    return passed, tested


def main():
    import itertools
    n = 4
    perms = Permutation.all_permutations(n)
    for i in range(1, n - 1):
        isobar = beta_isobaric_operator(i, beta)
        sref = (nh.one - (x[i] - x[i+1])*nh(simple_reflection(i)))
        divd = nh(simple_reflection(i))
        for perm1, perm2 in itertools.product(perms, repeat=2):
            f = Sx(perm1)
            g = Sx(perm2)
            lhs = isobar.apply(f * g)
            rhs = isobar.apply(f) * g + sref.apply(f) * (isobar.apply(g) + beta * g)
            diff =lhs - rhs
            if diff != Sx.zero:
                print("FAILED")
                print("i=", i, "f=", f, "g=", g)
                print("lhs=", lhs)
                print("rhs=", rhs)
                print("diff=", diff)


if __name__ == "__main__":
    main()
