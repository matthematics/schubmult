"""Quantum (single) Schubert polynomial multiplication.

Implements ``schubmult_q``/``schubmult_q_fast``: the product of a linear
combination of quantum Schubert polynomials ``S_u`` with a single ``S_v``,
returned as a coefficient dict ``{w: coeff}`` polynomial in the quantum
parameters ``q_1, q_2, ...``. Uses the same ``theta``/v-path recursion as
``schubmult.mult.single``, with the elementary-symmetric step generalized to
the quantum Pieri-type moves of ``elem_sym_perms_q`` (which may pick up a
factor of ``q`` when a Bruhat move is replaced by its quantum analogue).
``schubmult_q`` uses ``strict_theta`` (no repeated layer merging);
``schubmult_q_fast``/``_schubmult_q_fast_python`` uses ``medium_theta`` and
merges adjacent equal-length layers via ``double_elem_sym_q`` for speed.
"""

from functools import cached_property

import schubmult.utils.schub_lib as sss
from schubmult.combinatorics.permutation import (
    Permutation,
    uncode,
)
from schubmult.mult import _accel
from schubmult.symbolic import Add, Mul, Pow
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

logger = get_logger(__name__)


class _gvars:
    """Lazily-constructed default generating sets (avoids building symbols at import time)."""

    @cached_property
    def n(self):
        return 100

    @cached_property
    def var_x(self):
        return GeneratingSet("x")

    @cached_property
    def q_var(self):
        return GeneratingSet("q")


_vars = _gvars()


def single_variable(coeff_dict, varnum, var_q=_vars.q_var):
    """Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum`` (quantum Monk rule).

    Same structure as ``schubmult.mult.single.single_variable``, using
    ``elem_sym_perms_q`` so that some Bruhat moves carry a factor from ``var_q``.
    """
    ret = {}
    for u in coeff_dict:
        new_perms_k = sss.elem_sym_perms_q(u, 1, varnum, var_q)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = sss.elem_sym_perms_q(u, 1, varnum - 1, var_q)
        # print(f"{new_perms_k=}")
        for perm, udiff, mul_val in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u] * mul_val
        for perm, udiff, mul_val in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u] * mul_val
    return ret


def mult_poly_q(coeff_dict, poly, var_x=_vars.var_x, var_q=_vars.q_var):
    """Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

    Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
    single-variable leaves to ``single_variable``; mirrors ``mult_poly_py``.
    """
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    # logger.debug(f"{poly=} {type(poly)=} {list(var_x)}")
    # logger.debug(f"{[type(v) for v in var_x]}")
    if var_x.index(poly) != -1:
        # logger.debug(f"Found {var_x.index(poly)=}")
        # print("bang")
        return single_variable(coeff_dict, var_x.index(poly), var_q=var_q)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_q(ret, a, var_x, var_q=var_q)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_q(ret, base, var_x, var_q=var_q)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_q(coeff_dict, a, var_x, var_q=var_q))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def schubmult_q_fast(perm_dict, v, q_var=_vars.q_var):
    """Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

    Dispatches to the compiled ``schubmult_cpp`` kernel when available, falling
    back to ``_schubmult_q_fast_python`` (the ``medium_theta``-based recursion
    with merged equal-length layers) otherwise.

    Args:
        perm_dict: Mapping ``{Permutation: coeff}``.
        v: Permutation (or array-form list) indexing the quantum Schubert
            polynomial to multiply by.
        q_var: Generating set for the quantum parameters.

    Returns:
        dict: Coefficient dict ``{Permutation: coeff}``, polynomial in ``q_var``.
    """
    if _accel.available:
        ret = _accel.schubmult_q_fast(perm_dict, v, q_var)
        if ret is not None:
            return ret
    return _schubmult_q_fast_python(perm_dict, v, q_var)


def _schubmult_q_fast_python(perm_dict, v, q_var=_vars.q_var):
    """Pure-Python implementation of ``schubmult_q_fast``; see there for the contract."""
    if v.inv == 0:
        return perm_dict
    th = (~v).medium_theta()
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}

    thL = len(th)
    # if thL!=2 and len(set(thL))!=1:
    # raise ValueError("Not what I can do")
    vpathdicts = sss.compute_vpathdicts(th, vmu)
    # print(f"{vpathdicts=}")
    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            if index > 0 and th[index - 1] == th[index]:
                continue
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            if index < len(th) - 1 and th[index] == th[index + 1]:
                mx_th1 = 0
                for vp in vpathdicts[index + 1]:
                    for v2, vdiff, s in vpathdicts[index + 1][vp]:
                        mx_th1 = max(mx_th1, th[index + 1] - vdiff)
                newpathsums = {}
                for up in vpathsums:
                    newpathsums0 = {}
                    inv_up = up.inv
                    newperms = sss.double_elem_sym_q(up, mx_th, mx_th1, th[index], q_var)
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff2, s2 in vpathdicts[index][v]:
                            for up1, udiff1, mul_val1 in newperms:
                                if (up1, udiff1, mul_val1) not in newpathsums0:
                                    newpathsums0[(up1, udiff1, mul_val1)] = {}
                                if udiff1 + vdiff2 == th[index]:
                                    newpathsums0[(up1, udiff1, mul_val1)][v2] = (
                                        newpathsums0[(up1, udiff1, mul_val1)].get(
                                            v2,
                                            0,
                                        )
                                        + s2 * sumval * mul_val1
                                    )

                    for up1, udiff1, mul_val1 in newpathsums0:
                        for v in vpathdicts[index + 1]:
                            sumval = newpathsums0[(up1, udiff1, mul_val1)].get(v, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff2, s2 in vpathdicts[index + 1][v]:
                                for up2, udiff2, mul_val2 in newperms[(up1, udiff1, mul_val1)]:
                                    if up2 not in newpathsums:
                                        newpathsums[up2] = {}
                                    if udiff2 + vdiff2 == th[index + 1]:
                                        newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s2 * sumval * mul_val2
            else:
                newpathsums = {}
                for up in vpathsums:
                    inv_up = up.inv
                    newperms = sss.elem_sym_perms_q(
                        up,
                        min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                        th[index],
                        q_var,
                    )
                    for up2, udiff, mul_val in newperms:
                        if up2 not in newpathsums:
                            newpathsums[up2] = {}
                        for v in vpathdicts[index]:
                            sumval = vpathsums[up].get(v, 0)
                            if sumval == 0:
                                continue
                            for v2, vdiff, s in vpathdicts[index][v]:
                                if udiff + vdiff == th[index]:
                                    newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s * sumval * mul_val
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict(
            {ep: vpathsums[ep].get(toget, 0) for ep in vpathsums},
            ret_dict,
        )
    return ret_dict


def schubmult_q(perm_dict, v):
    """Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

    Reference (non-"fast") implementation: uses ``strict_theta`` and processes
    every layer individually (no merging of equal-length adjacent layers), so it
    is simpler but slower than ``schubmult_q_fast``. Results agree with
    ``schubmult_q_fast`` for all inputs.

    Args:
        perm_dict: Mapping ``{Permutation: coeff}``.
        v: Permutation (or array-form list) indexing the quantum Schubert
            polynomial to multiply by.

    Returns:
        dict: Coefficient dict ``{Permutation: coeff}``, polynomial in the
        default quantum parameters ``q``.
    """
    th = (~v).strict_theta()
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    if len(th) == 0:
        return perm_dict
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = sss.compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {u: {Permutation([]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = sss.elem_sym_perms_q(
                    up,
                    min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                    th[index],
                )
                for up2, udiff, mul_val in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            if udiff + vdiff == th[index]:
                                newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + s * sumval * mul_val
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict(
            {ep: vpathsums[ep].get(toget, 0) for ep in vpathsums},
            ret_dict,
        )
    return ret_dict
