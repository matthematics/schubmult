"""Multiplication kernels for (single) beta-Grothendieck polynomials.

* ``grothmult_py``: the product ``G_u * G_v`` in the ``G`` basis, by the same route as
  ``groth_double.grothmult_double`` specialized to ``y = z = 0`` -- expand ``G_v`` into
  Schubert polynomials via WC graphs, then push each ``S_{v'}`` through the v-path layers
  of ``theta(v'^{-1})`` with the closed-form K-Pieri coefficient ``groth_elem_sym_coeff``
  (a binomial times a power of ``beta``).
* ``single_variable_groth`` / ``mult_poly_groth``: multiplication by ``x_k`` (the
  non-equivariant K-theoretic Chevalley formula) and by arbitrary polynomials in ``x``.

The Chevalley coefficients are derived from M. Willems, "A Chevalley formula in
equivariant K-theory" (arXiv:math/0603220), Theorem 5 (the ordinary, non-equivariant
specialization of his equivariant Chevalley formula, Theorem 4). Willems indexes K-theory
classes O_w by the *dimension* of the Schubert variety, dual to the *codimension* indexing
used by Schubert/Grothendieck polynomials S_w/G_w; the w0-conjugation in
``_chevalley_ev_weights`` (``hat_w = w0*w`` going in, ``w0*v`` coming out) translates
between the two conventions. The beta-grading (beta^(d-1) per length difference
d = l(v)-l(w)) matches this codebase's beta-deformed Grothendieck polynomial
normalization (beta=0 recovers the classical Monk formula).
"""

from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import Add, Mul, Pow
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet_base
from schubmult.utils.perm_utils import add_perm_dict


def _chevalley_ev_weights(w, k, n):
    """Willems' ordinary K-theory Chevalley coefficients ``ev(v)`` for weight
    ``lambda = -e_k``, computed via the reduced word of ``w0(n) * w`` processed
    right-to-left through the Demazure (T^0) / reflection (T^1) operators of
    Willems' Theorem 2, then dualized back via w0. Returns ``{v: ev}``.
    """
    w0 = Permutation.w0(n)
    hat_w = w0 * w
    word = hat_w.code_word
    lam = [0] * n
    lam[k - 1] = -1
    state = {Permutation([]): {tuple(lam): 1}}
    for mu in reversed(word):
        s = Permutation.ref_product(mu)
        new_state = {}
        for perm, wdict in state.items():
            # epsilon = 1: reflection operator T^1_mu, perm picks up s on the left
            new_perm = s @ perm
            nd = new_state.setdefault(new_perm, {})
            for wv, c in wdict.items():
                wv2 = list(wv)
                wv2[mu - 1], wv2[mu] = wv2[mu], wv2[mu - 1]
                nd[tuple(wv2)] = nd.get(tuple(wv2), 0) + c
            # epsilon = 0: Demazure operator T^0_mu, perm unchanged
            nd0 = new_state.setdefault(perm, {})
            for wv, c in wdict.items():
                p = wv[mu - 1] - wv[mu]
                if p == 0:
                    continue
                if p > 0:
                    for jj in range(p):
                        wv2 = list(wv)
                        wv2[mu - 1] -= jj
                        wv2[mu] += jj
                        nd0[tuple(wv2)] = nd0.get(tuple(wv2), 0) + c
                else:
                    for jj in range(1, -p + 1):
                        wv2 = list(wv)
                        wv2[mu - 1] += jj
                        wv2[mu] -= jj
                        nd0[tuple(wv2)] = nd0.get(tuple(wv2), 0) - c
        state = new_state
    return {w0 * perm: ev for perm, wdict in state.items() if (ev := sum(wdict.values())) != 0}


def chevalley_x_k(w, k, beta, n=None):
    """Coefficients of ``x_k * G_w^(beta)`` in the Grothendieck basis, as a dict
    ``{v: coeff}`` (``w`` itself never appears: the self-term cancels identically).
    """
    w = Permutation(w)
    if n is None:
        n = max(len(w), k) + 1
    result = {}
    for v, ev in _chevalley_ev_weights(w, k, n).items():
        if v == w:
            continue
        d = v.inv - w.inv
        coeff = ((-1) ** d) * ev * beta ** (d - 1)
        if coeff != 0:
            result[v] = result.get(v, 0) + coeff
    return result


def single_variable_groth(coeff_dict, varnum, beta):
    """Multiply ``sum_u coeff_u G_u^(beta)`` by the single variable ``x_varnum``
    (Grothendieck Chevalley formula), via ``chevalley_x_k``.
    """
    ret = {}
    for u, coeff in coeff_dict.items():
        for v, c in chevalley_x_k(u, varnum, beta).items():
            ret[v] = ret.get(v, 0) + c * coeff
    return ret


def mult_poly_groth(coeff_dict, poly, var_x, beta):
    """Multiply ``sum_u coeff_u G_u^(beta)`` by an arbitrary polynomial ``poly`` in ``var_x``.

    Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
    single-variable leaves to ``single_variable_groth``.
    """
    if not isinstance(var_x, GeneratingSet_base):
        var_x = CustomGeneratingSet(var_x)
    if var_x.index(poly) != -1:
        return single_variable_groth(coeff_dict, var_x.index(poly), beta)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for a in poly.args:
            ret = mult_poly_groth(ret, a, var_x, beta)
        return ret
    if isinstance(poly, Pow):
        base = poly.args[0]
        exponent = int(poly.args[1])
        ret = coeff_dict
        for _ in range(exponent):
            ret = mult_poly_groth(ret, base, var_x, beta)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for a in poly.args:
            ret = add_perm_dict(ret, mult_poly_groth(coeff_dict, a, var_x, beta))
        return ret
    ret = {}
    for perm in coeff_dict:
        ret[perm] = poly * coeff_dict[perm]
    return ret


def groth_elem_sym_coeff(k, u1, u2, vdiff, beta):
    r"""Coefficient of ``G_{u2}`` in ``E_{k - vdiff, k}(x; 0) G_{u1} = e_{k - vdiff}(x_1..x_k) G_{u1}``.

    The ``y = z = 0`` specialization of ``groth_double._groth_elem_sym_frac``.  Sort the
    window positions ``j <= k`` by the fate of ``u1(j)`` in ``u2``: *fixed*, *left*
    (reappears at an earlier window position) or *out* (leaves the window or moves
    right); with ``F = #fixed``, ``L = #left``, ``m = L + #out`` movers and
    ``d = l(u2) - l(u1)`` the closed form

        beta^(d - m) (-beta)^L E_{n - q, n}( (-)y_fixed, (-1/beta)^L ; z ),   n = F + L, q = vdiff,

    collapses at ``y = z = 0`` (the fixed alphabet entries become ``0``) to

        beta^(d - m) (-beta)^(L - p) binom(L, p),   p = n - q,   0 <= p <= L,

    and to ``0`` otherwise.  At ``beta = 0`` only ``d = m``, ``p = L`` survive, i.e. the
    classical Pieri rule ``e_p(x_1..x_k) S_u = sum_{u ->_k w} S_w`` with ``p = k - d``
    (there ``L = 0`` for chains with distinct lower indices, so ``p = F = k - d``).
    """
    from math import comb

    d = u2.inv - u1.inv
    window2 = [u2[j] for j in range(k)]
    fixed = left = 0
    for j in range(k):
        value = u1[j]
        if window2[j] == value:
            fixed += 1
        elif value in window2 and window2.index(value) < j:
            left += 1
    movers = k - fixed
    if d < movers:
        return 0
    p = fixed + left - vdiff
    if p < 0 or p > left:
        return 0
    return beta ** (d - movers) * (-beta) ** (left - p) * comb(left, p)


def _groth_schub_vpath_mul(perm_dict, v, beta):
    """``sum_u coeff_u G_u(x) * S_v(x)`` in the ``G`` basis: the ``y = z = 0`` case of
    ``groth_double._groth_schub_vpath_mul``, with the same v-path layers and the same
    marked-chain supports, but integer/``beta``-polynomial coefficients throughout.
    """
    from schubmult.combinatorics.permutation import uncode
    from schubmult.mult.groth_double import _top_block_support
    from schubmult.utils.schub_lib import compute_vpathdicts

    v = Permutation(v)
    th = list((~v).theta())
    while th and th[-1] == 0:
        th.pop()
    if not th:
        return dict(perm_dict)
    vmu = v * uncode(th)
    vpathdicts = compute_vpathdicts(tuple(th), vmu)
    ret_dict = {}
    for u, val in perm_dict.items():
        u = Permutation(u)
        vpathsums = {u: {Permutation([1, 2]): val}}
        for index in range(len(th)):
            k = th[index]
            newpathsums = {}
            for up, sums in vpathsums.items():
                for up2 in _top_block_support(up, k) | {up}:
                    for v_iter, steps in vpathdicts[index].items():
                        sumval = sums.get(v_iter)
                        if sumval is None or sumval == 0:
                            continue
                        for v2, vdiff, s in steps:
                            coeff = groth_elem_sym_coeff(k, up, up2, vdiff, beta)
                            if coeff == 0:
                                continue
                            bucket = newpathsums.setdefault(up2, {})
                            bucket[v2] = bucket.get(v2, 0) + s * sumval * coeff
            vpathsums = newpathsums
        for ep, sums in vpathsums.items():
            value = sums.get(vmu)
            if value is not None and value != 0:
                ret_dict[ep] = ret_dict.get(ep, 0) + value
    return {w: c for w, c in ret_dict.items() if c != 0}


def grothmult_py(perm_dict, v, beta=None):
    r"""Multiply (single) Grothendieck polynomials, mirroring ``schubmult_py``.

    Computes the expansion of ``sum_u coeff_u G_u(x) * G_v(x)`` in the basis ``{G_w(x)}``
    and returns it as ``{w: coeff_w}`` with coefficients polynomial in ``beta``.

    Same method as ``groth_double.grothmult_double`` specialized to ``y = z = 0``: expand
    ``G_v`` into Schubert polynomials (``groth_elem_as_schub_dict``, via WC graphs) and
    push each ``S_{v'}`` through the v-path kernel ``_groth_schub_vpath_mul``, whose
    per-layer coefficients are the binomial closed form ``groth_elem_sym_coeff``.
    """
    from schubmult.abc import beta as default_beta
    from schubmult.symbolic.poly.schub_poly import groth_elem_as_schub_dict

    if beta is None:
        beta = default_beta
    v = Permutation(v)
    perm_dict = {Permutation(key): value for key, value in perm_dict.items()}
    if v.inv == 0:
        return perm_dict
    ret = {}
    for vprime, coeff in groth_elem_as_schub_dict(v, beta).items():
        for w, value in _groth_schub_vpath_mul(perm_dict, vprime, beta).items():
            ret[w] = ret.get(w, 0) + coeff * value
    from schubmult.symbolic import expand

    return {w: expand(c) for w, c in ret.items() if expand(c) != 0}
