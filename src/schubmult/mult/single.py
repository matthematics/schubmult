"""Ordinary (single) Schubert polynomial multiplication.

Implements the ``schubmult_py`` kernel: the product of a linear combination of
Schubert polynomials ``S_u`` (given as a dict ``{u: coeff}``) with a single
Schubert polynomial ``S_v``, returned as a coefficient dict ``{w: coeff}``.

The algorithm is the recursive "v-path" / transition method used throughout
``schubmult``: write ``theta = (~v).theta()`` for the dominant weakly-decreasing
vector bounding ``v``'s Lehmer code, let ``mu = uncode(theta)`` and
``vmu = v * mu``, and process the entries of ``theta`` one at a time, tracking
Bruhat-chain "v-paths" from ``vmu`` down to the identity (``compute_vpathdicts``)
alongside chains of elementary-symmetric moves on the ``u`` side
(``elem_sym_perms``). Accumulating consistent pairs of chains and reading off
the coefficient landing on ``vmu`` gives the product.
"""

from functools import cached_property

from schubmult.combinatorics.permutation import (
    Permutation,
    uncode,
)
from schubmult.mult import _accel
from schubmult.symbolic import Add, Mul, Pow
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
from schubmult.utils.logging import get_logger, init_logging
from schubmult.utils.perm_utils import (
    add_perm_dict,
)
from schubmult.utils.schub_lib import (
    compute_vpathdicts,
    elem_sym_perms,
    elem_sym_perms_op,
)

init_logging(debug=False)
logger = get_logger(__name__)

class _gvars:
    """Lazily-constructed default generating set (avoids building symbols at import time)."""

    @cached_property
    def var_x(self):
        return GeneratingSet("x")


_vars = _gvars()


def single_variable(coeff_dict, varnum):
    """Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum``.

    Uses the classical Monk rule: ``x_k * S_u = sum S_{u t_{ij}}`` over Bruhat
    covers ``u t_{ij}`` with ``i <= k < j`` (added) minus those with ``j <= k < i``
    (subtracted), via ``elem_sym_perms(u, 1, varnum)``.

    Args:
        coeff_dict: Mapping ``{Permutation: coeff}``.
        varnum: 1-indexed variable index ``k``.

    Returns:
        dict: The updated coefficient dict ``{Permutation: coeff}``.
    """
    ret = {}
    for u in coeff_dict:
        new_perms_k = elem_sym_perms(u, 1, varnum)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = elem_sym_perms(u, 1, varnum - 1)
        for perm, udiff in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u]
        for perm, udiff in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u]
    return ret


# TODO: if need indexes, CustomGeneratingSet
def mult_poly_py(coeff_dict, poly, var_x=_vars.var_x):
    """Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

    Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``; each single
    variable leaf is dispatched to ``single_variable``, and any other leaf just
    scales every coefficient.

    Args:
        coeff_dict: Mapping ``{Permutation: coeff}``.
        poly: Symbolic polynomial expression in the variables of ``var_x``.
        var_x: Generating set identifying the ``x`` variables (default ``x``).

    Returns:
        dict: The updated coefficient dict ``{Permutation: coeff}``.
    """
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)

    if var_x.index(poly) != -1:
        # print(f"{poly=} {var_x._symbols_arr=} {var_x._symbols_arr.index(poly)=}")
        return single_variable(coeff_dict, var_x.index(poly))

    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_py(ret, a, var_x)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_py(ret, base, var_x)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_py(coeff_dict, a, var_x))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def schubmult_py(perm_dict, v):
    """Multiply ``sum_u coeff_u S_u(x)`` by the (ordinary) Schubert polynomial ``S_v``.

    Dispatches to the compiled ``schubmult_cpp`` kernel when available and the
    permutations fit within its ``MAXN``, falling back to the pure-Python
    implementation otherwise.

    Args:
        perm_dict: Mapping ``{Permutation: coeff}`` (integer coefficients).
        v: Permutation (or array-form list) indexing the Schubert polynomial to
            multiply by.

    Returns:
        dict: Coefficient dict ``{Permutation: coeff}`` for the product.
    """
    if _accel.available:
        ret = _accel.schubmult_py(perm_dict, v)
        if ret is not None:
            return ret
    return _schubmult_py_python(perm_dict, v)


def _schubmult_py_python(perm_dict, v):
    """Pure-Python implementation of ``schubmult_py``; see there for the contract."""
    v = Permutation(v)
    # print(f"{v=}")
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    vpathdicts = compute_vpathdicts(th, vmu)
    mx_th = [0 for i in range(len(th))]
    for index in range(len(th)):
        for vp in vpathdicts[index]:
            for v2, vdiff, s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {Permutation(u): {Permutation([1, 2]): val}}

        for index in range(len(th)):
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms(
                    up,
                    min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                    th[index],
                )
                # print(f"{up=}")
                for vp in vpathsums[up]:
                    # print(f"{vp=} {type(vp)=} {hash(vp)=}")
                    # print(f"{vpathsums[up]=} {vpathdicts[index]=}")
                    sumval = vpathsums[up][vp]
                    if sumval == 0:
                        continue
                    for v2, vdiff, s in vpathdicts[index][vp]:
                        addsumval = s * sumval
                        for up2, udiff in newperms:
                            if vdiff + udiff == th[index]:
                                if up2 not in newpathsums:
                                    newpathsums[up2] = {}
                                newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + addsumval
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict

def schubmult_py_down(perm_dict, v):
    """Divided-difference ("down") variant of ``_schubmult_py_python``.

    Same v-path recursion but built from ``elem_sym_perms_op`` (Bruhat *descents*)
    instead of ``elem_sym_perms``, used for the down/dual side of the transition
    recursion rather than ordinary multiplication.

    Args:
        perm_dict: Mapping ``{Permutation: coeff}``.
        v: Permutation to multiply by.

    Returns:
        dict: Coefficient dict ``{Permutation: coeff}``.
    """
    v = Permutation(v)
    # print(f"{v=}")
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    vpathdicts = compute_vpathdicts(th, vmu)
    mx_th = [0 for i in range(len(th))]
    for index in range(len(th)):
        for vp in vpathdicts[index]:
            for v2, vdiff, s in vpathdicts[index][vp]:
                mx_th[index] = max(mx_th[index], th[index] - vdiff)

    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {Permutation(u): {Permutation([1, 2]): val}}

        for index in range(len(th)):
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms_op(
                    up,
                    min(mx_th[index], inv_mu - inv_vmu - (inv_up - inv_u)),
                    th[index],
                )
                # print(f"{up=}")
                for vp in vpathsums[up]:
                    # print(f"{vp=} {type(vp)=} {hash(vp)=}")
                    # print(f"{vpathsums[up]=} {vpathdicts[index]=}")
                    sumval = vpathsums[up][vp]
                    if sumval == 0:
                        continue
                    for v2, vdiff, s in vpathdicts[index][vp]:
                        addsumval = s * sumval
                        for up2, udiff in newperms:
                            if vdiff + udiff == th[index]:
                                if up2 not in newpathsums:
                                    newpathsums[up2] = {}
                                newpathsums[up2][v2] = newpathsums[up2].get(v2, 0) + addsumval
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schub_coprod_py(perm, indices):
    """Coproduct of ``S_perm`` restricted to the variable split named by ``indices``.

    Computes the expansion of the (single) Schubert polynomial coproduct
    ``Delta_{indices}(S_perm) = sum (firstperm, secondperm) -> coeff`` by
    multiplying the Grassmannian permutation for ``indices`` against ``perm``
    via ``schubmult_py`` and splitting each resulting permutation's window into
    its first ``N`` and remaining ``len(perm) - N`` values.

    Args:
        perm: Permutation (or array-form list) to take the coproduct of.
        indices: Iterable of 1-indexed positions selecting the variable split.

    Returns:
        dict: Mapping ``{(firstperm, secondperm): coeff}``.
    """
    mperm = perm
    indices = sorted(indices)
    ret_dict = {}
    k = len(indices)
    n = len(mperm)
    kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    max_required = max([kcd[i] + i for i in range(len(kcd))])
    kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
    N = len(kcd)
    kperm = ~(uncode(kcd2))
    coeff_dict = {kperm: 1}
    #logger.debug(f"{kperm.code}*{mperm.code}")
    coeff_dict = schubmult_py(coeff_dict, mperm)

    inv_kperm = kperm.inv
    inverse_kperm = ~kperm
    # total_sum = 0
    for perm, val in coeff_dict.items():
        if val == 0:
            continue
        # pperm = [*perm]
        downperm = perm * inverse_kperm
        if downperm.inv == perm.inv - inv_kperm:
            flag = True
            for i in range(N):
                if downperm[i] > N:
                    flag = False
                    break
            if not flag:
                continue
            firstperm = Permutation(downperm[0:N])
            secondperm = Permutation([downperm[i] - N for i in range(N, len(downperm))])
            ret_dict[(firstperm, secondperm)] = val
    return ret_dict
