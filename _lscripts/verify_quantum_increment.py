"""Pin down the quantum increment in the induction step of the Pieri proof.

The single-variable quantum Monk (Lemma 3.2) splits as
    (x_k - z_j) S_v = CLASSICAL[ (y_{v(k)}-z_j)S_v + sum_raise sign S_{v t_{kp}} ]
                    + DROP_UP  [ sum_{b>k, drop} +q_{kb} S_{v t_{kb}} ]
                    + DROP_DOWN[ sum_{a<k, drop} -q_{ak} S_{v t_{ak}} ]

The elementary recursion is
    S_u E^q_{p,k} = (x_k - z_{k-p+1}) S_u E^q_{p-1,k-1} + S_u E^q_{p,k-1}
                    + q_{k-1} S_u E^q_{p-2,k-2}

We test two structural claims that would make the write-up clean:

  (Q1) DROP_DOWN( S_u E^q_{p-1,k-1} )  +  q_{k-1} S_u E^q_{p-2,k-2}  ==  0
       (the invalid a<k drop terms are cancelled exactly by the recursion's quantum term)

  (Q2) CLASSICAL( S_u E^q_{p-1,k-1} ) + DROP_UP( S_u E^q_{p-1,k-1} ) + S_u E^q_{p,k-1}
       ==  S_u E^q_{p,k}
       (equivalently: everything except DROP_DOWN, plus E_{p,k-1}, already gives the target)

If BOTH hold for all (u,p,k), the proof splits cleanly into a classical skeleton
(CLASSICAL + E_{p,k-1} -> length-RAISE ->_k terms, ports from arXiv:2401.11060) plus a
self-contained quantum piece (DROP_UP -> length-DROP ->_k terms carrying q_{kb}; and the
DROP_DOWN/recursion cancellation Q1).

Run:  conda activate schubmult_312 && python _lscripts/verify_quantum_increment.py 4
"""

import sys

from schubmult import *  # noqa: F401,F403
from schubmult.abc import x
from schubmult.combinatorics.permutation import Permutation
from schubmult.rings.schubert import QDSx
from schubmult.symbolic import S, expand, prod
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.symbolic.symmetric_polynomials import QFactorialElemSym

z_gs = GeneratingSet("z")
q_gs = GeneratingSet("q")


def q_ab(a, b):
    return prod([q_gs[s] for s in range(a, b)])


def product_Eq(Su, p, k):
    ring = Su.ring
    if p < 0 or p > k:
        return ring.zero
    if p == 0:
        return ring.from_dict(Su)
    Epk = QFactorialElemSym(p, k, [x[i] for i in range(1, k + 1)], [z_gs[i] for i in range(1, k - p + 2)])
    return ring.elem_mul(Su, Epk)


def monk_split(elem, i, j, N):
    """Return (classical, drop_up, drop_down) pieces of (x_i - z_j)*elem via Lemma 3.2.

    classical  = diagonal + all length-RAISE reflection terms
    drop_up    = length-DROP terms with reflection partner b>i (carry +q_{i,b})
    drop_down  = length-DROP terms with reflection partner a<i (carry -q_{a,i})
    """
    ring = elem.ring
    yg = ring.coeff_genset
    classical = ring.zero
    drop_up = ring.zero
    drop_down = ring.zero
    for v, c in elem.items():
        classical += (c * (yg[v[i - 1]] - z_gs[j])) * ring(v)
        for pp in range(1, N + 1):
            if pp == i:
                continue
            vt = v.swap(i - 1, pp - 1)
            dl = vt.inv - v.inv
            sign = 1 if pp > i else -1
            if dl == 1:
                classical += (c * sign) * ring(vt)
            elif dl == -2 * abs(pp - i) + 1:
                qf = q_ab(min(i, pp), max(i, pp))
                if pp > i:
                    drop_up += (c * sign * qf) * ring(vt)
                else:
                    drop_down += (c * sign * qf) * ring(vt)
    return classical, drop_up, drop_down


def is_zero(elem):
    return all(expand(coeff) == 0 for coeff in elem.values())


def elems_equal(a, b):
    return is_zero(a - b)


def run(n):
    N = 2 * n
    perms = list(Permutation.all_permutations(n))
    total = 0
    q1_fail = []
    q2_fail = []
    for u in perms:
        Su = QDSx(u)
        ring = Su.ring
        for k in range(1, n):
            for p in range(1, k + 1):
                j = k - p + 1
                classical, drop_up, drop_down = monk_split(product_Eq(Su, p - 1, k - 1), k, j, N)
                # (Q1): drop_down + q_{k-1} E_{p-2,k-2} == 0
                q1 = drop_down + q_gs[k - 1] * product_Eq(Su, p - 2, k - 2)
                # (Q2): classical + drop_up + E_{p,k-1} == E_{p,k}
                q2 = classical + drop_up + product_Eq(Su, p, k - 1)
                target = product_Eq(Su, p, k)
                total += 1
                if not is_zero(q1):
                    q1_fail.append((tuple(u), p, k))
                if not elems_equal(q2, target):
                    q2_fail.append((tuple(u), p, k))
    print(f"n={n}: tested {total} steps;  Q1(dropdown+recursion==0) failures={len(q1_fail)};  Q2(rest==target) failures={len(q2_fail)}")
    for m in q1_fail[:10]:
        print("  Q1 FAIL:", m)
    for m in q2_fail[:10]:
        print("  Q2 FAIL:", m)
    return len(q1_fail) == 0 and len(q2_fail) == 0


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    ok = run(n)
    sys.exit(0 if ok else 1)
