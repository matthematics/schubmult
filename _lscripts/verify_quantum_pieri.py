"""Numerically verify the conjectured quantum-double Pieri formula (Theorem 1.2 of
"A Molev-Sagan type formula for quantum double Schubert polynomials").

Statement:
    S^q_u(x;y) * E^q_{p,k}(x;z)
        = sum_{u ->_k w} q_k(u,w) * E_{p-n_k, k-n_k}(y_{P_k(u,w)}; z) * S^q_w(x;y)

We test it as:  (honest generic quantum-double product of two double Schuberts in
DIFFERENT coefficient sets y and z)  ==  (the ring's `elem_mul` Pieri shortcut).

    LHS = S^q_u(x;y) * S^q_{c_{p,k}}(x;z)         [ generic pair product ]
    RHS = elem_mul( S^q_u(x;y), E^q_{p,k}(x;z) )  [ conjectured Pieri expansion ]

where c_{p,k} = s_{k-p+1} s_{k-p+2} ... s_k is the cycle with S_{c_{p,k}} = E_{p,k}.

Run:  conda activate schubmult_312 && python _lscripts/verify_quantum_pieri.py 3
"""

import sys

from schubmult import *  # noqa: F401,F403  (initialize package to avoid circular import in lazy submodules)
from schubmult.abc import x, z
from schubmult.combinatorics.permutation import Permutation
from schubmult.rings.schubert import QDSx
from schubmult.symbolic import expand
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import QFactorialElemSym

z_gs = GeneratingSet("z")


def cycle_cpk(p, k):
    """Return the permutation c_{p,k} = s_{k-p+1} s_{k-p+2} ... s_k as a Permutation."""
    return Permutation.ref_product(*range(k - p + 1, k + 1))


def elems_equal(a, b):
    """Compare two ring elements robustly."""
    if hasattr(a, "almosteq"):
        try:
            if a.almosteq(b):
                return True
        except Exception:  # noqa: BLE001
            pass
    diff = a - b
    for v in diff.values():
        if expand(v) != 0:
            return False
    return True


def run(n):
    perms = list(Permutation.all_permutations(n))
    total = 0
    mismatches = []
    for u in perms:
        Su = QDSx(u)
        ring_y = Su.ring
        for k in range(1, n):
            for p in range(1, k + 1):
                c = cycle_cpk(p, k)
                # honest product of two quantum double Schuberts in different coeff sets
                lhs = Su * QDSx(c, z_gs)
                # Pieri shortcut: multiply by E^q_{p,k}(x; z)
                Epk = QFactorialElemSym(p, k, [x[i] for i in range(1, k + 1)], [z[i] for i in range(1, k + 2 - p)])
                rhs = ring_y.elem_mul(Su, Epk)
                total += 1
                if not elems_equal(lhs, rhs):
                    mismatches.append((tuple(u), p, k))
    print(f"n={n}: tested {total} (u,p,k) triples; mismatches={len(mismatches)}")
    for m in mismatches[:20]:
        print("  MISMATCH:", m)
    return len(mismatches) == 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    ok = run(n)
    sys.exit(0 if ok else 1)
