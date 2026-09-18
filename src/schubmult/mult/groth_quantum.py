r"""Multiplication kernel for (single) quantum beta-Grothendieck polynomials.

``grothmult_q`` computes ``G^q_u(x) * G^q_v(x)`` in the ``G^q`` basis: the ``y = z = 0``
specialization of ``groth_quantum_double.grothmult_q_double``, by the same route as
``groth.grothmult_py`` -- expand ``G_v`` into Schubert polynomials, push each ``S_{v'}``
through the strict-theta v-path layers, and use the quantum K-Pieri support
``quantum_pieri_chains`` with the binomial closed form ``groth_elem_sym_coeff`` evaluated at
the quantum chain length and weighted by ``q^D``.

``G^q_v = Q(G_v)`` is the Lenart--Maeno quantization (``groth_quantum_double.lm_quantize``);
at ``beta = -1`` these are the quantum Grothendieck polynomials of Lenart--Maeno with
``Q_j = q_j``, and the ``e_p`` Pieri rule ``grothmult_q_pieri`` is then equivalent to the
Naito--Sagaki quantum K Pieri theorem (arXiv:2211.01578): the sign ``(-1)^{len - p}`` and the
marking count ``#Mark`` there are the ``beta^{len - m} (-beta)^{L - p'} binom(L, p')`` here,
summed over the ``e_p <-> G_{c[k,p]}`` triangular change of basis.  Conjectural in general;
see ``groth_quantum_double`` for the evidence.
"""

from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.mult.groth import groth_elem_sym_coeff
from schubmult.mult.groth_quantum_double import _qmon, quantum_pieri_chains
from schubmult.symbolic import expand
from schubmult.symbolic.poly.schub_poly import _vars
from schubmult.utils.perm_utils import add_perm_dict
from schubmult.utils.schub_lib import compute_vpathdicts

__all__ = ["grothmult_q", "grothmult_q_dict", "grothmult_q_pieri"]


def grothmult_q_pieri(coeff_dict, p, k, beta=None, q_var=None):
    r"""Multiply ``sum_u coeff_u G^q_u(x)`` by ``Q(e_p(x_1..x_k))``.

    Coefficient of ``G^q_w``: ``q^D * groth_elem_sym_coeff(k, u, w, k - p, beta, length)``
    over the quantum K-Pieri support; ``p = k`` is the non-equivariant quantum top block.
    """
    from schubmult.abc import beta as default_beta

    if beta is None:
        beta = default_beta
    if q_var is None:
        q_var = _vars.q_var
    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        if p == 0:
            ret[u] = ret.get(u, 0) + val
            continue
        for w, (length, dvec) in quantum_pieri_chains(u, k).items():
            coeff = groth_elem_sym_coeff(k, u, w, k - p, beta, length=length)
            if coeff != 0:
                ret[w] = ret.get(w, 0) + val * _qmon(dvec, q_var) * coeff
    return ret


def _qgroth_schub_vpath_mul(perm_dict, v, beta, q_var):
    """``sum_u coeff_u G^q_u(x) * Q(S_v(x))`` in the ``G^q`` basis.

    ``groth._groth_schub_vpath_mul`` with ``strict_theta`` (products along a strict-theta
    v-path are standard elementary monomials, on which the quantization is multiplicative),
    ``quantum_pieri_chains`` as the support, and the chain length in the coefficient.
    """
    v = Permutation(v)
    th = list((~v).strict_theta())
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
        for index, k in enumerate(th):
            newpathsums = {}
            for up, sums in vpathsums.items():
                for up2, (length, dvec) in quantum_pieri_chains(up, k).items():
                    qmon = _qmon(dvec, q_var)
                    for v_iter, steps in vpathdicts[index].items():
                        sumval = sums.get(v_iter)
                        if sumval is None or sumval == 0:
                            continue
                        for v2, vdiff, s in steps:
                            coeff = groth_elem_sym_coeff(k, up, up2, vdiff, beta, length=length)
                            if coeff == 0:
                                continue
                            bucket = newpathsums.setdefault(up2, {})
                            bucket[v2] = bucket.get(v2, 0) + s * sumval * qmon * coeff
            vpathsums = newpathsums
        for ep, sums in vpathsums.items():
            value = sums.get(vmu)
            if value is not None and value != 0:
                ret_dict[ep] = ret_dict.get(ep, 0) + value
    return {w: c for w, c in ret_dict.items() if c != 0}


def grothmult_q(perm_dict, v, beta=None, q_var=None):
    r"""Multiply (single) quantum Grothendieck polynomials, mirroring ``schubmult_q``.

    Returns the expansion of ``sum_u coeff_u G^q_u(x) * G^q_v(x)`` in the basis ``{G^q_w(x)}``
    as ``{w: coeff_w}``, polynomial in ``beta`` and ``q``.  ``G_v`` is expanded into Schubert
    polynomials by ``groth_elem_as_schub_dict`` (quantization is linear over ``beta``) and
    each ``S_{v'}`` pushed through ``_qgroth_schub_vpath_mul``.  ``q = 0`` is ``grothmult_py``;
    ``beta = 0`` is ``schubmult_q``.
    """
    from schubmult.abc import beta as default_beta
    from schubmult.symbolic.poly.schub_poly import groth_elem_as_schub_dict

    if beta is None:
        beta = default_beta
    if q_var is None:
        q_var = _vars.q_var
    v = Permutation(v)
    perm_dict = {Permutation(key): value for key, value in perm_dict.items()}
    if v.inv == 0:
        return perm_dict
    ret = {}
    for vprime, coeff in groth_elem_as_schub_dict(v, beta).items():
        for w, value in _qgroth_schub_vpath_mul(perm_dict, vprime, beta, q_var).items():
            ret[w] = ret.get(w, 0) + coeff * value
    out = {w: expand(c) for w, c in ret.items()}
    return {w: c for w, c in out.items() if c != 0}


def grothmult_q_dict(perm_dict1, perm_dict2, beta=None, q_var=None):
    """Product of two coefficient dicts: ``sum_v coeff2_v grothmult_q(perm_dict1, v, ...)``."""
    ret = {}
    for v, coeff in perm_dict2.items():
        ret = add_perm_dict(ret, {w: coeff * value for w, value in grothmult_q(perm_dict1, v, beta, q_var).items()})
    return ret
