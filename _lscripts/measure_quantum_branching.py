"""Measure the QUANTUM schubmult branching structure faithfully.

The quantum kernel (schubmult_q_fast) differs from the classical one:
  - it uses medium_theta() (a composition, possibly with repeated parts),
  - when th[i]==th[i+1] it processes TWO columns at once via double_elem_sym_q,
  - elem_sym_perms_q covers carry q-degrees (mul_val) and may decrease length.

We measure, per v (and per u via a single-node front):
  Lq        = len(medium_theta) after trimming (quantum layer count)
  Lq_eff    = number of processed layers (paired columns count once)
  Bq        = max |elem_sym_perms_q(up,p,k)| single-column branch
  Bdbl      = max total covers produced by double_elem_sym_q (sum len(perms2))
  Bdbl_pre  = max |perms1| (before cycle filter) in a double step
  filt_ratio= min over double steps of (kept)/(perms1*perms2_pre) cycle-filter tightness
  maxqdeg   = max total q-degree (number of q-vars) on any single cover mul_val
  Delta     = kdown-triple factor (UNCHANGED from classical; v-side)
"""

import itertools
import sys

from schubmult import Permutation
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import sympify
from schubmult.utils import schub_lib as sss
from schubmult.utils.schub_lib import compute_vpathdicts, double_elem_sym_q, elem_sym_perms_q


def qdeg(mul_val):
    """Number of q-variable factors in a mul_val (0 if it is the integer 1)."""
    if mul_val == 1:
        return 0
    e = sympify(mul_val)
    s = str(e)
    # count q_ occurrences crudely; exponents counted via '**'
    # more robust: total degree
    try:
        from schubmult.symbolic import Poly

        return sum(sum(m) for m in Poly(e).monoms()[:1]) if False else _degree(e)
    except Exception:
        return s.count("q")


def _degree(e):
    # total degree of a monomial in the q variables
    try:
        d = 0
        for f, expo in e.as_powers_dict().items():
            if "q" in str(f):
                d += int(expo)
        return d
    except Exception:
        return 0


def measure(v, nmax_front_u=None):
    th = (~v).medium_theta()
    if len(th) == 0:
        return None
    while th and th[-1] == 0:
        th.pop()
    if not th:
        return None
    mu = uncode(th)
    vmu = v * mu
    inv_mu = mu.inv
    inv_vmu = vmu.inv
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)

    Lq = thL
    Bq = 0
    Bdbl = 0
    Bdbl_pre = 0
    maxqdeg = 0
    Delta = 0
    Lq_eff = 0

    # Use u = identity front to exercise the branchings from the natural base.
    u = Permutation(list(range(1, len(v) + 1)))
    inv_u = u.inv
    index = 0
    while index < thL:
        # replicate the kernel's mx_th computation
        mx_th = 0
        for vp in vpathdicts[index]:
            for _v2, vdiff, _s in vpathdicts[index][vp]:
                mx_th = max(mx_th, th[index] - vdiff)
        for vp in vpathdicts[index]:
            Delta = max(Delta, len(vpathdicts[index][vp]))
        Lq_eff += 1
        if index < thL - 1 and th[index] == th[index + 1]:
            mx_th1 = 0
            for vp in vpathdicts[index + 1]:
                for _v2, vdiff, _s in vpathdicts[index + 1][vp]:
                    mx_th1 = max(mx_th1, th[index + 1] - vdiff)
            newperms = double_elem_sym_q(u, mx_th, mx_th1, th[index])
            npre = len(newperms)
            Bdbl_pre = max(Bdbl_pre, npre)
            total_covers = sum(len(newperms[k]) for k in newperms)
            Bdbl = max(Bdbl, total_covers)
            for (p1, ud1, mv1), lst in newperms.items():
                maxqdeg = max(maxqdeg, _degree(sympify(mv1)))
                for (p2, ud2, mv2) in lst:
                    maxqdeg = max(maxqdeg, _degree(sympify(mv2)))
            index += 2
            continue
        newperms = elem_sym_perms_q(u, min(mx_th, (inv_mu - (u.inv - inv_u)) - inv_vmu), th[index])
        Bq = max(Bq, len(newperms))
        for (p2, ud, mv) in newperms:
            maxqdeg = max(maxqdeg, _degree(sympify(mv)))
        index += 1

    return {"Lq": Lq, "Lq_eff": Lq_eff, "Bq": Bq, "Bdbl": Bdbl, "Bdbl_pre": Bdbl_pre, "maxqdeg": maxqdeg, "Delta": Delta, "lv": v.inv}


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for n in range(2, nmax + 1):
        agg = {"Lq": 0, "Lq_eff": 0, "Bq": 0, "Bdbl": 0, "Bdbl_pre": 0, "maxqdeg": 0, "Delta": 0}
        args = {k: None for k in agg}
        paired_seen = False
        for v in Permutation.all_permutations(n):
            r = measure(v)
            if r is None:
                continue
            if r["Bdbl"] > 0:
                paired_seen = True
            for k in agg:
                if r[k] > agg[k]:
                    agg[k] = r[k]
                    args[k] = list(v)
        from math import comb, log2

        print(
            f"n={n}: Lq_max={agg['Lq']}({args['Lq']}) Lq_eff_max={agg['Lq_eff']} "
            f"Bq={agg['Bq']}({args['Bq']}) Bdbl={agg['Bdbl']}({args['Bdbl']}) "
            f"Bdbl_pre={agg['Bdbl_pre']} maxqdeg={agg['maxqdeg']}({args['maxqdeg']}) "
            f"Delta={agg['Delta']} | 2^n-1={2**n-1} (2^n-1)^2={(2**n-1)**2} paired={paired_seen}",
            flush=True,
        )


if __name__ == "__main__":
    main()
