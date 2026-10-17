"""Prototype: schubmult_double with a mod-p shadow that prunes states whose partial sum is (probabilistically) zero."""

import random
import sys
import time

from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import sympify
from schubmult.symbolic.common_polys import elem_sym_func
from schubmult.symbolic.poly.variables import GeneratingSet
from schubmult.utils.perm_utils import add_perm_dict
from schubmult.utils.schub_lib import compute_vpathdicts, elem_sym_perms

P = (1 << 61) - 1
zero = sympify(0)


def make_shadow(seed):
    rnd = random.Random(seed)
    point = {}
    memo = {}

    def ev(e):
        v = memo.get(e)
        if v is not None:
            return v
        if e.is_Integer:
            v = int(e) % P
        elif e.is_Symbol:
            v = point.setdefault(e, rnd.randrange(1, P))
        elif e.is_Add:
            v = sum(ev(a) for a in e.args) % P
        elif e.is_Mul:
            v = 1
            for a in e.args:
                v = v * ev(a) % P
        elif e.is_Pow:
            b, x = e.args
            v = pow(ev(b), int(x), P)
        elif e.is_Rational:
            v = int(e.p) * pow(int(e.q), -1, P) % P
        else:
            raise TypeError(type(e))
        memo[e] = v
        return v

    return ev


def schubmult_double_shadow(perm_dict, v, var2, var3, prune, stats, verify=None):
    """Same DP as _schubmult_double_python; each state carries (Expr, value mod p).

    With ``prune`` states whose value is 0 are dropped at the end of each level; with ``verify`` (a
    list) the dropped expressions are collected for an exact check.
    """
    ev = make_shadow(1)
    perm_dict = {Permutation(k): vv for k, vv in perm_dict.items()}
    v = Permutation(v)
    th = (~v).theta()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu, inv_mu = vmu.inv, mu.inv
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {u: {Permutation([1, 2]): (val, ev(val))}}
        for index in range(thL):
            mx_th = max(th[index] - vdiff for vp in vpathdicts[index] for _, vdiff, _ in vpathdicts[index][vp])
            pending = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms(up, min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu), th[index])
                for up2, udiff in newperms:
                    tgt = pending.setdefault(up2, {})
                    for v_iter, (sumval, sumv) in vpathsums[up].items():
                        for v2, vdiff, s in vpathdicts[index][v_iter]:
                            esf = elem_sym_func(th[index], index + 1, up, up2, v_iter, v2, udiff, vdiff, var2, var3)
                            if esf == 0:
                                continue
                            term = s * sumval * esf
                            termv = (s * sumv * ev(esf)) % P
                            tgt.setdefault(v2, []).append((term, termv))
            newpathsums = {}
            n_states = n_zero = 0
            for up2, d in pending.items():
                for v2, terms in d.items():
                    total = sum((t for t, _ in terms), zero)
                    totalv = sum(tv for _, tv in terms) % P
                    if total == 0:
                        continue
                    n_states += 1
                    if totalv == 0:
                        n_zero += 1
                        if prune:
                            if verify is not None:
                                verify.append(total)
                            continue
                    newpathsums.setdefault(up2, {})[v2] = (total, totalv)
            stats.append((index, n_states, n_zero))
            vpathsums = newpathsums
        ret_dict = add_perm_dict({ep: vpathsums[ep][vmu][0] for ep in vpathsums if vmu in vpathsums[ep]}, ret_dict)
    return ret_dict


if __name__ == "__main__":
    y, z = GeneratingSet("y"), GeneratingSet("z")
    u, v = Permutation([4, 1, 6, 5, 2, 3]), Permutation([8, 1, 7, 6, 2, 3, 5, 4])
    same = len(sys.argv) > 1 and sys.argv[1] == "same"
    v3 = y if same else z
    for prune in (False, True):
        stats = []
        pruned = [] if prune else None
        t = time.perf_counter()
        out = schubmult_double_shadow({u: sympify(1)}, v, y, v3, prune, stats, pruned)
        dt = time.perf_counter() - t
        print(f"prune={prune}: {dt:.2f}s, {len(out)} output terms")
        for index, n, nz in stats:
            print(f"  level {index}: {n} live states, {nz} of them zero mod p")
        if pruned:
            from schubmult.symbolic import expand

            t = time.perf_counter()
            bad = sum(1 for e in pruned if expand(e) != 0)
            print(f"  exact verification of {len(pruned)} pruned states by expansion: {time.perf_counter() - t:.2f}s, {bad} false zeros")
