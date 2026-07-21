"""Verify that the natural induction on k CLOSES for the quantum-double Pieri formula.

This is the decisive test for the proof of Theorem 1.2 of
"A Molev-Sagan type formula for quantum double Schubert polynomials".

The proof is an induction on the level k, driven by the single-variable recursion for
the quantum factorial elementary polynomial (schub_poly.py `elem_sym_poly_q`):

    E^q_{p,k}(x;z) = (x_k - z_{k-p+1}) E^q_{p-1,k-1}(x;z)
                     + E^q_{p,k-1}(x;z)
                     + q_{k-1} E^q_{p-2,k-2}(x;z)                              (*)

Assuming the Pieri expansions of  S^q_u * E^q_{p-1,k-1},  S^q_u * E^q_{p,k-1},  and
S^q_u * E^q_{p-2,k-2}  are known (the induction hypothesis at levels k-1, k-2), we must
multiply the FIRST of these by (x_k - z_{k-p+1}) and add.  The ONLY tool used for that
multiplication is the single-variable quantum Monk formula (Lemma 3.2), applied with i=k:

    (x_k - z_j) S^q_v(x;y) = (y_{v(k)} - z_j) S^q_v(x;y)
        + sum_{p != k, l(v t_{kp}) = l(v)+1}          sign(p-k)        S^q_{v t_{kp}}(x;y)
        + sum_{p != k, l(v t_{kp}) = l(v)-2|p-k|+1}   sign(p-k) q_{kp} S^q_{v t_{kp}}(x;y)

We implement Lemma 3.2 from scratch (`monk_mul`) and the recursion (*), then check that

    monk_mul( S_u*E^q_{p-1,k-1}, i=k, j=k-p+1 )
      + S_u*E^q_{p,k-1}
      + q_{k-1} * S_u*E^q_{p-2,k-2}
    ==  S_u * E^q_{p,k}

for all u in S_n and all 1 <= p <= k < n.  The degree-(k-1)/(k-2) expansions and the
degree-k target are all produced by the ring's (trusted, already-verified) `elem_mul`.

A pass means: GIVEN the quantum Monk formula, the induction closes exactly -- there is
no leftover term and no missing term.  This is precisely what the written proof must show.

Run:  conda activate schubmult_312 && python _lscripts/verify_pieri_induction.py 4
"""

import sys

from schubmult import *  # noqa: F401,F403  (init package to avoid circular import)
from schubmult.abc import x
from schubmult.combinatorics.permutation import Permutation
from schubmult.rings.schubert import QDSx
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import QFactorialElemSym

z_gs = GeneratingSet("z")
q_gs = GeneratingSet("q")


def q_ab(a, b):
    """q_{a,b} = q_a q_{a+1} ... q_{b-1}  (a < b)."""
    return prod([q_gs[s] for s in range(a, b)])


def product_Eq(Su, p, k):
    """Return the ring element  S^q_u * E^q_{p,k}(x; z_1..z_{k-p+1})  via trusted elem_mul.

    Handles the base pieces of recursion (*):  E^q_{0,k}=1, and E^q_{p,k}=0 for p<0 or p>k.
    """
    ring = Su.ring
    if p < 0 or p > k:
        return ring.zero
    if p == 0:
        return ring.from_dict(Su)  # E^q_{0,k} = 1
    Epk = QFactorialElemSym(p, k, [x[i] for i in range(1, k + 1)], [z_gs[i] for i in range(1, k - p + 2)])
    return ring.elem_mul(Su, Epk)


def monk_mul(elem, i, j, N):
    """Multiply a ring element by (x_i - z_j) using the single-variable quantum Monk (Lemma 3.2).

    p ranges over positions 1..N (p != i); N must exceed the support of every basis perm.
    """
    ring = elem.ring
    yg = ring.coeff_genset
    out = ring.zero
    for v, c in elem.items():
        # diagonal term (y_{v(i)} - z_j) S^q_v
        out += (c * (yg[v[i - 1]] - z_gs[j])) * ring(v)
        for pp in range(1, N + 1):
            if pp == i:
                continue
            vt = v.swap(i - 1, pp - 1)
            dl = vt.inv - v.inv
            sign = 1 if pp > i else -1
            if dl == 1:
                out += (c * sign) * ring(vt)
            elif dl == -2 * abs(pp - i) + 1:
                out += (c * sign * q_ab(min(i, pp), max(i, pp))) * ring(vt)
    return out


def elems_equal(a, b):
    diff = a - b
    for coeff in diff.values():
        if expand(coeff) != 0:
            return False
    return True


def run(n):
    N = 2 * n  # generous support bound for reflections t_{k,b}
    perms = list(Permutation.all_permutations(n))
    total = 0
    mismatches = []
    for u in perms:
        Su = QDSx(u)
        for k in range(1, n):
            for p in range(1, k + 1):
                j = k - p + 1  # z-index used by recursion (*): varl2[k-p] = z_{k-p+1}
                # induction combination (RHS of (*), pushed through S^q_u)
                comb = monk_mul(product_Eq(Su, p - 1, k - 1), k, j, N)
                comb += product_Eq(Su, p, k - 1)
                comb += q_gs[k - 1] * product_Eq(Su, p - 2, k - 2)
                # target: level-k expansion
                target = product_Eq(Su, p, k)
                total += 1
                if not elems_equal(comb, target):
                    mismatches.append((tuple(u), p, k))
    print(f"n={n}: tested {total} (u,p,k) induction steps; mismatches={len(mismatches)}")
    for m in mismatches[:20]:
        print("  MISMATCH:", m)
    return len(mismatches) == 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    ok = run(n)
    sys.exit(0 if ok else 1)
