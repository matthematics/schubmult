"""beta-Grothendieck Chevalley formula: multiplication of a (single) Grothendieck
polynomial by a bare x_k variable, i.e. x_k * G_w^(beta).

Non-equivariant (y=0) for now; ``GrothendieckRing`` has no coefficient/y genset yet.

Derived from M. Willems, "A Chevalley formula in equivariant K-theory"
(arXiv:math/0603220), Theorem 5 (the ordinary, non-equivariant specialization of
his equivariant Chevalley formula, Theorem 4). Willems indexes K-theory classes
O_w by the *dimension* of the Schubert variety, dual to the *codimension*
indexing used by Schubert/Grothendieck polynomials S_w/G_w; the w0-conjugation
below (``hat_w = w0*w`` going in, ``w0*v`` coming out) translates between the two
conventions. The beta-grading (beta^(d-1) per length difference d = l(v)-l(w))
matches this codebase's beta-deformed Grothendieck polynomial normalization
(beta=0 recovers the classical double Schubert Monk formula). Calibrated against
grothendieck_poly()/to_groth() (see session notes).
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
    ret = {}
    for u, coeff in coeff_dict.items():
        for v, c in chevalley_x_k(u, varnum, beta).items():
            ret[v] = ret.get(v, 0) + c * coeff
    return ret


def mult_poly_groth(coeff_dict, poly, var_x, beta):
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
