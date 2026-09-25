"""Double Schubert polynomial multiplication.

Implements the ``schubmult_double`` kernel: the product of a linear
combination of double Schubert polynomials ``S_u(x, var2)`` with a single
``S_v(x, var3)``, returned as a coefficient dict ``{w: coeff}`` of polynomials
in ``var2``/``var3``. Uses the same ``theta``/``vmu``/v-path recursion as
``schubmult.mult.single`` (see that module), replacing the plain elementary
symmetric contribution with ``elem_sym_func``, which carries the secondary
variables ``var2``/``var3``.

Also provides the "alt"/"from_elems" variants (building the product one
descent-pulled variable at a time via ``pull_out_var``, generic to any choice
of elementary-symmetric-like function), ``nilhecke_mult`` (nilHecke ring
multiplication), and ``schub_coprod_double`` (the double Schubert coproduct).
"""

from bisect import bisect_left
from functools import cache

from schubmult.combinatorics.permutation import (
    Permutation,
    uncode,
)
from schubmult.mult import _accel
from schubmult.symbolic import Add, Mul, Pow, S, expand, expand_func, sympify
from schubmult.symbolic.common_polys import _vars, efficient_subs, elem_func_func_mul, elem_sym_func, elem_sym_poly
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict
from schubmult.utils.schub_lib import (
    compute_vpathdicts,
    elem_sym_perms,
    elem_sym_perms_op,
    elem_sym_positional_perms,
    pull_out_var,
)

zero = sympify(0)

logger = get_logger(__name__)


def count_sorted(mn, tp):
    """Count occurrences of ``tp`` in the sorted sequence ``mn`` via binary search."""
    # BUG: unused elsewhere in the codebase; the while loop never advances `index`, so it
    # infinite-loops whenever mn[index] == tp. Left as-is/uncovered pending a fix decision.
    index = bisect_left(mn, tp)
    ct = 0
    if mn[index] == tp:
        while index < len(mn) and mn[index] == tp:
            ct += 1
    return ct


# def E(p, k, varl=None):
#     return elem_sym_poly(p, k, _vars.var1[1:], varl)


def single_variable(coeff_dict, varnum, var2=None):
    """Multiply ``sum_u coeff_u S_u(x, var2)`` by the single variable ``x_varnum``.

    Equivariant Monk rule: the diagonal term contributes ``var2[u(varnum)]``
    (localization of ``x_varnum`` at ``u``) and the off-diagonal terms are the
    same Bruhat-cover moves as the ordinary (non-equivariant) ``single_variable``
    in ``schubmult.mult.single``.

    Args:
        coeff_dict: Mapping ``{Permutation: coeff}``.
        varnum: 1-indexed variable index ``k``.
        var2: Secondary (``y``) generating set.

    Returns:
        dict: The updated coefficient dict.
    """
    ret = {}
    for u in coeff_dict:
        if varnum - 1 < len(u):
            ret[u] = ret.get(u, 0) + var2[u[varnum - 1]] * coeff_dict[u]
        else:
            ret[u] = ret.get(u, 0) + var2[varnum] * coeff_dict[u]
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


def single_variable_down(coeff_dict, varnum, var2=None):
    """Down (descent) variant of ``single_variable``, using ``elem_sym_perms_op``."""
    ret = {}
    for u in coeff_dict:
        if varnum - 1 < len(u):
            ret[u] = ret.get(u, 0) + var2[u[varnum - 1]] * coeff_dict[u]
        else:
            ret[u] = ret.get(u, 0) + var2[varnum] * coeff_dict[u]
        new_perms_k = elem_sym_perms_op(u, 1, varnum)
        new_perms_km1 = []
        if varnum > 1:
            new_perms_km1 = elem_sym_perms_op(u, 1, varnum - 1)
        for perm, udiff in new_perms_k:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) + coeff_dict[u]
        for perm, udiff in new_perms_km1:
            if udiff == 1:
                ret[perm] = ret.get(perm, 0) - coeff_dict[u]
    return ret


def mult_poly_double(coeff_dict, poly, var_x=None, var_y=None):
    """Multiply ``sum_u coeff_u S_u(x, var_y)`` by an arbitrary polynomial ``poly`` in ``var_x``.

    Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
    single-variable leaves to ``single_variable``; mirrors ``mult_poly_py`` with
    the extra secondary alphabet ``var_y``.
    """
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    if var_x.index(poly) != -1:
        return single_variable(coeff_dict, var_x.index(poly), var_y)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_double(ret, a, var_x, var_y)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_double(ret, base, var_x, var_y)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_double(coeff_dict, a, var_x, var_y))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def mult_poly_double_alt(coeff_dict, poly, var_x=None, var_y=None):
    """Variant of ``mult_poly_double`` that folds each factor via ``schubmult_double_dict``
    instead of ``single_variable``, so ``poly`` is only ever expanded one variable/factor
    at a time in the ``S_v`` basis rather than left as a raw scalar multiplier.
    """
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    if var_x.index(poly) != -1:
        return single_variable(coeff_dict, var_x.index(poly), var_y)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            s_d = mult_poly_double_alt({Permutation([]): S.One}, a, var_x, var_y)
            ret = schubmult_double_dict(ret, s_d, var_y, var_y)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        s_d = mult_poly_double_alt({Permutation([]): S.One}, base, var_x, var_y)
        for i in range(int(exponent)):
            ret = schubmult_double_dict(ret, s_d, var_y, var_y)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_double_alt(coeff_dict, a, var_x, var_y))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


# def mult_poly_symy(coeff_dict, poly, var_x=_vars.sympy_var1, var_y=_vars.sympy_var2):


def mult_poly_down(coeff_dict, poly):
    """Down (descent) variant of ``mult_poly_double``, using ``single_variable_down``
    and the fixed default alphabet ``_vars.var1``.
    """
    # BUG: unused elsewhere in the codebase; single_variable_down is called below without a
    # var2, so it crashes (None is not subscriptable) on any polynomial containing a real
    # _vars.var1 symbol. Left as-is/uncovered pending a fix decision.
    if poly in _vars.var1:
        return single_variable_down(coeff_dict, _vars.var1.index(poly))
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_down(ret, a)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for i in range(int(exponent)):
            ret = mult_poly_down(ret, base)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_down(coeff_dict, a))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def nilhecke_mult(coeff_dict1, coeff_dict2):
    """NilHecke ring product of ``coeff_dict1`` (polynomial coefficients) and
    ``coeff_dict2`` (permutation coefficients acting as divided-difference operators).

    For each ``w`` in ``coeff_dict2`` its coefficient polynomial is pushed through
    ``mult_poly_down`` against ``coeff_dict1``, and each resulting permutation ``v``
    is right-multiplied by ``w`` whenever that multiplication is length-additive.

    Returns:
        dict: Coefficient dict ``{Permutation: coeff}``.
    """
    # BUG: unused elsewhere in the codebase; `v1 = [*v]` unpacks the Permutation into a plain
    # list, so `v1 * w1` below fails (can't multiply a list by a Permutation) whenever did_mul
    # is nonempty. Left as-is/uncovered pending a fix decision.
    ret = {}
    for w in coeff_dict2:
        w1 = w
        inv_w1 = w1.inv
        poly = coeff_dict2[w]
        did_mul = mult_poly_down(coeff_dict1, poly)
        for v in did_mul:
            v1 = [*v]
            addperm = v1 * w1
            if addperm.inv == v1.inv + inv_w1:
                toadd = addperm
                ret[toadd] = ret.get(toadd, 0) + did_mul[v]
    return ret


@cache
def schubmult_double_pair(perm1, perm2, var2=None, var3=None):
    """``schubmult_double`` specialized to a single ``perm1`` with coefficient 1, cached."""
    return schubmult_double({perm1: S.One}, perm2, var2, var3)


@cache
def schubmult_double_pair_generic(perm1, perm2):
    """``schubmult_double_pair`` with the fixed generic secondary alphabets ``_vars.var_g1``/``_vars.var_g2``."""
    return schubmult_double({perm1: S.One}, perm2, _vars.var_g1, _vars.var_g2)


@cache
def schubmult_double_pair_generic_alt(perm1, perm2):
    """Like ``schubmult_double_pair_generic`` but computed via ``schubmult_double_alt_from_elems``
    with the factorial elementary symmetric function, then expanded/simplified.
    """
    from schubmult.symbolic.symmetric_polynomials import FactorialElemSym

    return {k: expand_func(expand(v)) for k, v in schubmult_double_alt_from_elems({perm1: S.One}, perm2, _vars.var_g1, _vars.var_g2, elem_func=FactorialElemSym).items()}


def schubmult_double_dict(perm_dict1, perm_dict2, var2=None, var3=None):
    """Multiply two coefficient dicts of double Schubert polynomials together.

    Computes ``(sum_u coeff1_u S_u(x, var2)) * (sum_v coeff2_v S_v(x, var3))``
    by summing ``schubmult_double(perm_dict1, v, var2, var3)`` scaled by
    ``coeff2_v`` over ``v`` in ``perm_dict2``.
    """
    ret = {}
    for k, v in perm_dict2.items():
        ret = add_perm_dict(ret, {k2: v2 * v for k2, v2 in schubmult_double(perm_dict1, k, var2, var3).items()})
    return ret


def schubmult_double(perm_dict, v, var2=None, var3=None):
    """Multiply ``sum_u coeff_u S_u(x, var2)`` by the double Schubert polynomial ``S_v(x, var3)``.

    Dispatches to the compiled ``schubmult_cpp`` kernel when available (and both
    secondary alphabets are given), falling back to the pure-Python
    implementation ``_schubmult_double_python`` otherwise.

    Args:
        perm_dict: Mapping ``{Permutation: coeff}``.
        v: Permutation (or array-form list) indexing the Schubert polynomial to
            multiply by.
        var2: Secondary alphabet attached to ``perm_dict``'s permutations.
        var3: Secondary alphabet attached to ``v``.

    Returns:
        dict: Coefficient dict ``{Permutation: coeff}`` (polynomials in ``var2``/``var3``).
    """
    if _accel.available and var2 is not None and var3 is not None:
        ret = _accel.schubmult_double(perm_dict, v, var2, var3)
        if ret is not None:
            return ret
    return _schubmult_double_python(perm_dict, v, var2, var3)


def _schubmult_double_python(perm_dict, v, var2=None, var3=None):
    """Pure-Python implementation of ``schubmult_double``; see there for the contract."""
    perm_dict = {Permutation(k): vv for k, vv in perm_dict.items()}
    v = Permutation(v)
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0:
        return perm_dict
    if th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms(
                    up,
                    min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                    th[index],
                )
                for up2, udiff in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v_iter in vpathdicts[index]:
                        sumval = vpathsums[up].get(v_iter, zero)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v_iter]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                zero,
                            ) + s * sumval * elem_sym_func(
                                th[index],
                                index + 1,
                                up,
                                up2,
                                v_iter,
                                v2,
                                udiff,
                                vdiff,
                                var2,
                                var3,
                            )
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({Permutation(ep): vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schubmult_double_alt(perm_dict, v, var2=None, var3=None, index=1):
    """Alternate double Schubert product, built by peeling one variable of ``~v`` at a
    time via ``pull_out_var`` instead of the ``theta``/v-path recursion.

    Multiplies ``sum_u coeff_u S_u(x, var2)`` by ``S_v(x, var3)``, recursing on
    ``~new_v`` with the elementary symmetric factor coming from
    ``elem_sym_positional_perms`` at each step.
    """
    if v.inv == 0:
        return perm_dict
        # ret = S.Zero
    ret_dict = {}
    L = pull_out_var(1, ~v)
    for index_list, new_v in L:
        interim_dict = {}
        for u, val in perm_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_sym_poly(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        ret_dict = add_perm_dict(ret_dict, schubmult_double_alt(interim_dict, ~new_v, var2, var3, index + 1))
    return ret_dict


# forwards backwards
def schubmult_double_alt_from_elems_forwards(perm_dict, v, var2=None, var3=None, index=1, elem_func=None):
    """``schubmult_double_alt`` generalized to an arbitrary elementary-symmetric-like
    ``elem_func(p, k, x_vars, y_vars)``, processing variables of ``~v`` from the first
    pulled-out index forward.
    """
    if v.inv == 0:
        return perm_dict
        # ret = S.Zero
    ret_dict = {}
    L = pull_out_var(1, ~v)
    for index_list, new_v in L:
        interim_dict = {}
        for u, val in perm_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        ret_dict = add_perm_dict(ret_dict, schubmult_double_alt_from_elems_forwards(interim_dict, ~new_v, var2, var3, index + 1, elem_func))
    return ret_dict


# backwards mul after
# def schubmult_double_alt_from_elems(perm_dict, v, var2=None, var3=None, elem_func=None):
#     if v.inv == 0:
#         return perm_dict
#     ret_dict = {}
#     index = max((~v).descents()) + 1
#     L = pull_out_var(index, ~v)
#     for index_list, new_v in L:
#         interim_dict = {}
#         for u, val in perm_dict.items():
#             new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
#             for new_perm, p, sgn in new_perms:
#                 interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
#                     len(index_list) - p,
#                     len(index_list) - p,
#                     [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
#                     [var3[index]],
#                 )
#         ret_dict = add_perm_dict(ret_dict, schubmult_double_alt_from_elems(interim_dict, ~new_v, var2, var3, elem_func))
#     return ret_dict


# backwards mul before
def schubmult_double_alt_from_elems_backwards(perm_dict, v, var2=None, var3=None, elem_func=None):
    """Like ``schubmult_double_alt_from_elems_forwards`` but processing ``~v``'s pulled-out
    variables from the last descent backward, multiplying the elementary-symmetric
    factor in *before* recursing (dispatches to the compiled kernel when available).
    """
    if _accel.available and var2 is not None and var3 is not None and elem_func is not None:
        ret = _accel.schubmult_double_alt_from_elems(perm_dict, v, var2, var3, elem_func)
        if ret is not None:
            return ret
    return _schubmult_double_alt_from_elems_backwards_python(perm_dict, v, var2, var3, elem_func)


def _schubmult_double_alt_from_elems_backwards_python(perm_dict, v, var2=None, var3=None, elem_func=None):
    """Pure-Python implementation of ``schubmult_double_alt_from_elems_backwards``."""
    if v.inv == 0:
        return perm_dict
    ret_dict = {}
    index = max((~v).descents()) + 1
    L = pull_out_var(index, ~v)
    _cache = {}
    for index_list, new_v in L:
        if new_v not in _cache:
            _cache[new_v] = schubmult_double_alt_from_elems_backwards(perm_dict, ~new_v, var2, var3, elem_func)
        start_dict = _cache[new_v]
        # start_dict = perm_dict
        interim_dict = {}
        for u, val in start_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        ret_dict = add_perm_dict(ret_dict, interim_dict)
        # ret_dict = add_perm_dict(ret_dict,schubmult_double_alt_from_elems_backwards(interim_dict, ~new_v, var2, var3, elem_func))
    return ret_dict


def schubmult_double_alt_from_elems_backwards_backwards(perm_dict, v, var2=None, var3=None, elem_func=None):
    """Variant of ``_schubmult_double_alt_from_elems_backwards_python`` without the
    per-``new_v`` memoization cache, recursing on the interim dict instead of the
    original ``perm_dict`` at each pulled-out variable.
    """
    if v.inv == 0:
        return perm_dict
    ret_dict = {}
    index = max((~v).descents()) + 1
    L = pull_out_var(index, ~v)
    # _cache = {}
    for index_list, new_v in L:
        # if new_v not in _cache:
        #     _cache[new_v] = schubmult_double_alt_from_elems_backwards(perm_dict, ~new_v, var2, var3, elem_func)
        start_dict = perm_dict
        # start_dict = perm_dict
        interim_dict = {}
        for u, val in start_dict.items():
            new_perms = elem_sym_positional_perms(u, len(index_list), *index_list)
            for new_perm, p, sgn in new_perms:
                interim_dict[new_perm] = interim_dict.get(new_perm, S.Zero) + sgn * val * elem_func(
                    len(index_list) - p,
                    len(index_list) - p,
                    [var2[new_perm[i - 1]] for i in index_list if new_perm[i - 1] == u[i - 1]],
                    [var3[index]],
                )
        # ret_dict = add_perm_dict(ret_dict, interim_dict)
        ret_dict = add_perm_dict(ret_dict, schubmult_double_alt_from_elems_backwards_backwards(interim_dict, ~new_v, var2, var3, elem_func))
    return ret_dict


schubmult_double_alt_from_elems = schubmult_double_alt_from_elems_backwards


def schubmult_double_from_elems(perm_dict, v, var2=None, var3=None, elem_func=None):
    """``schubmult_double`` generalized to an arbitrary elementary-symmetric-like
    ``elem_func``, via the ``theta``/v-path recursion (rather than ``pull_out_var``).

    Dispatches to the compiled kernel when available, falling back to
    ``_schubmult_double_from_elems_python``.
    """
    if _accel.available and var2 is not None and var3 is not None and elem_func is not None:
        ret = _accel.schubmult_double_from_elems(perm_dict, v, var2, var3, elem_func)
        if ret is not None:
            return ret
    return _schubmult_double_from_elems_python(perm_dict, v, var2, var3, elem_func)


def _schubmult_double_from_elems_python(perm_dict, v, var2=None, var3=None, elem_func=None):
    """Pure-Python implementation of ``schubmult_double_from_elems``."""
    perm_dict = {Permutation(k): v for k, v in perm_dict.items()}
    v = Permutation(v)
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0:
        return perm_dict
    if th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    inv_vmu = vmu.inv
    inv_mu = mu.inv
    ret_dict = {}
    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        inv_u = u.inv
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                inv_up = up.inv
                newperms = elem_sym_perms(
                    up,
                    min(mx_th, (inv_mu - (inv_up - inv_u)) - inv_vmu),
                    th[index],
                )
                for up2, udiff in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, 0)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                0,
                            ) + s * sumval * elem_func_func_mul(
                                th[index],
                                index + 1,
                                up,
                                up2,
                                v,
                                v2,
                                udiff,
                                vdiff,
                                var2,
                                var3,
                                elem_func=elem_func,
                            )
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({Permutation(ep): vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schubmult_double_down(perm_dict, v, var2=None, var3=None):
    """Down (descent) variant of ``_schubmult_double_python``, using ``elem_sym_perms_op``."""
    vn1 = ~v
    th = vn1.theta()
    if len(th) == 0 or th[0] == 0:
        return perm_dict
    mu = uncode(th)
    vmu = v * mu
    ret_dict = {}

    while th[-1] == 0:
        th.pop()
    thL = len(th)
    vpathdicts = compute_vpathdicts(th, vmu)
    for u, val in perm_dict.items():
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(thL):
            mx_th = 0
            for vp in vpathdicts[index]:
                for v2, vdiff, s in vpathdicts[index][vp]:
                    mx_th = max(mx_th, th[index] - vdiff)
            newpathsums = {}
            for up in vpathsums:
                newperms = elem_sym_perms_op(up, mx_th, th[index])
                for up2, udiff in newperms:
                    if up2 not in newpathsums:
                        newpathsums[up2] = {}
                    for v in vpathdicts[index]:
                        sumval = vpathsums[up].get(v, zero)
                        if sumval == 0:
                            continue
                        for v2, vdiff, s in vpathdicts[index][v]:
                            newpathsums[up2][v2] = newpathsums[up2].get(
                                v2,
                                zero,
                            ) + s * sumval * elem_sym_func(
                                th[index],
                                index + 1,
                                up2,
                                up,
                                v,
                                v2,
                                udiff,
                                vdiff,
                                var2,
                                var3,
                            )
            vpathsums = newpathsums
        toget = vmu
        ret_dict = add_perm_dict({ep: vpathsums[ep].get(toget, 0) for ep in vpathsums}, ret_dict)
    return ret_dict


def schub_coprod_double(mperm, indices, var2=None, var3=None):
    """Coproduct of the double Schubert polynomial ``S_mperm`` restricted to the
    variable split named by ``indices``.

    Analogue of ``schub_coprod_py``: multiplies the Grassmannian permutation for
    ``indices`` against ``mperm`` (via ``schubmult_double`` with a merged ``2N``
    variable alphabet), splits each resulting permutation's window, and
    substitutes the merged alphabet back to ``var2``/``var3``.

    Args:
        mperm: Permutation (or array-form list) to take the coproduct of.
        indices: Iterable of 1-indexed positions selecting the variable split.
        var2: Secondary alphabet for the first factor's variables.
        var3: Secondary alphabet for the second factor's variables.

    Returns:
        dict: Mapping ``{(firstperm, secondperm): coeff}``.
    """
    indices = sorted(indices)
    subs_dict_coprod = {}
    k = len(indices)
    n = len(mperm)
    kcd = [indices[i] - i - 1 for i in range(len(indices))] + [n + 1 - k for i in range(k, n)]
    max_required = max([kcd[i] + i for i in range(len(kcd))])
    kcd2 = kcd + [0 for i in range(len(kcd), max_required)] + [0]
    N = len(kcd)
    kperm = ~uncode(kcd2)
    inv_kperm = kperm.inv
    vn = GeneratingSet("soible")

    for i in range(1, N * 2 + 1):
        if i <= N:
            subs_dict_coprod[vn[i]] = var2[i]
        else:
            subs_dict_coprod[vn[i]] = var3[i - N]

    coeff_dict = {kperm: 1}
    coeff_dict = schubmult_double(coeff_dict, mperm, vn, var2)

    inverse_kperm = ~kperm

    ret_dict = {}
    for perm in coeff_dict:
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

            val = efficient_subs(sympify(coeff_dict[perm]), subs_dict_coprod)

            key = (firstperm, secondperm)
            ret_dict[key] = val

    return ret_dict
