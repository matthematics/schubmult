r"""K-theoretic Monk formula for double Grothendieck polynomials.

Implements Lenart--Postnikov, *Affine Weyl groups in K-theory and representation
theory* (arXiv:math/0309207), Corollary 8.2 (the :math:`K_T`-Monk formula), in
type :math:`A` and transported to the ``schubmult`` conventions.

Reconciling the conventions
---------------------------
``schubmult`` has no exponentials: it uses the multiplicative formal group law
``a (+) b = a + b + beta*a*b`` with formal inverse ``(-)b = -b/(1 + beta*b)``.
The double Grothendieck polynomial of a simple reflection is

.. math::

    \mathfrak{G}_{s_k}(x, y)
        = \frac{1}{\beta}\Bigl(\prod_{i=1}^{k}(1 + \beta x_i)(1 + \beta y_i) - 1\Bigr),

which at ``beta = -1`` is precisely the class ``1 - x^{w_0(omega_k)} e^{-omega_k}``
of Lemma 8.1(a).  The dictionary between the paper and this package is

* torus characters: ``x^{eps_i}  <->  1 + beta*y_i`` (so that ``(x^{eps_a - eps_b} - 1)/beta``
  is the root ``y_a (-) y_b``, matching ``DoubleGrothendieckRing.exp_root``);
* basis elements: ``[O_{X_{w_0 w}}]  <->  (-beta)^{l(w)} G_w``, since the paper
  indexes structure sheaves by the *dimension* of the Schubert variety whereas
  ``G_w`` has lowest term ``S_w`` of degree ``l(w)`` (codimension indexing).

Under ``u -> w_0 u`` the saturated *decreasing* chains of Corollary 8.2 become
saturated *increasing* chains, the reflections ``t_{ij}`` are unchanged, and the
character prefactor ``x^{nu(J)} = x^{w_0(omega_k) - u(omega_k)}`` (constant in
``J`` because ``omega_k`` is minuscule in type ``A``) becomes
``prod_{i<=k} (1 + beta*y_i)/(1 + beta*y_{u(i)})``.  Rescaling by
``(-beta)^{l(.)}`` turns the signs ``(-1)^{|J|}`` into powers ``beta^{|J|}`` and
yields

.. math::

    \mathfrak{G}_u(x, y)\,\mathfrak{G}_{s_k}(x, z)
        = \frac{1}{\beta}\Bigl(\Theta_u \sum_J \beta^{|J|}\,
          \mathfrak{G}_{u\,r_J}(x, y) - \mathfrak{G}_u(x, y)\Bigr),
    \qquad
    \Theta_u = \prod_{i=1}^{k}\frac{1 + \beta z_i}{1 + \beta y_{u(i)}},

the sum being over the subsets ``J`` of a reduced ``(-omega_k)``-chain of
reflections whose reflections build a saturated increasing Bruhat chain from
``u`` (the empty subset included).  The two secondary alphabets are handled by
``G_{s_k}(x, z) = C G_{s_k}(x, y) + (C - 1)/beta`` with
``C = prod_{i<=k}(1 + beta*z_i)/(1 + beta*y_i)``, which is exactly what turns
the ``y_i`` of ``x^{nu}`` into the ``z_i`` of ``Theta_u``.

By Corollary 15.4 of the same paper a reduced ``(-omega_k)``-chain of
reflections in type ``A_{n-1}`` is

    ``t_{1,n}, t_{1,n-1}, ..., t_{1,k+1}, t_{2,n}, ..., t_{2,k+1}, ..., t_{k,k+1}``.
"""

from schubmult.abc import beta as _default_beta
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import Add, Mul, Pow, S, sympify, sympify_sympy
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet_base
from schubmult.utils.perm_utils import add_perm_dict

__all__ = [
    "epsilon_chain",
    "grothmult_double",
    "monk_chain",
    "mult_poly_groth_double",
    "single_variable_groth",
]


def monk_chain(k, n):
    """Reduced ``(-omega_k)``-chain of reflections in ``A_{n-1}`` (Cor. 15.4).

    Returned as the ordered tuple of transpositions ``(i, j)``, ``i <= k < j``.
    """
    return tuple((i, j) for i in range(1, k + 1) for j in range(n, k, -1))


def epsilon_chain(k, n):
    """Reduced ``(-eps_k)``-chain of reflections in ``A_{n-1}`` (Cor. 15.4).

    ``(-omega_k)``-chains only see the roots ``alpha_{ij}`` with ``i <= k < j``,
    which is why ``monk_chain`` multiplies by the whole product
    ``prod_{i<=k}(1 + beta*x_i)``.  Isolating a single ``x_k`` needs the weight
    ``eps_k = omega_k - omega_{k-1}`` instead, whose chain also involves the
    roots ``alpha_{ik}`` with ``i < k``.

    Each entry is ``(a, b, level)`` for the affine reflection ``s_{alpha_ab, level}``:

    * ``(i, k, 0)`` for ``i = k-1, ..., 1`` -- positive roots ``alpha_{ik}``;
    * ``(k, j, 1)`` for ``j = n, ..., k+1`` -- negative roots ``alpha_{jk}``, which
      are the ones contributing a sign ``(-1)^{n(J)}`` in Theorem 6.1.
    """
    return tuple((i, k, 0) for i in range(k - 1, 0, -1)) + tuple((k, j, 1) for j in range(n, k, -1))


def _one_plus_beta_x_terms(u, k, var2, beta, n):
    r"""Expand ``(1 + beta*x_k) G_u(x, var2)`` as ``{w: coeff}``.

    ``1 + beta*x_k`` is the class ``e^{-eps_k}``, so this is Theorem 6.1 for
    ``lambda = -eps_k``.  Writing ``J`` for a subset of ``epsilon_chain(k, n)``
    whose reflections form a saturated increasing Bruhat chain
    ``u = w_0 < w_1 < ... < w_s = w``, the transported coefficient is

        (-1)^{n(J)} (-beta)^{|J|} / (1 + beta*var2[w(k)])
            * prod_{level-1 steps} (1 + beta*var2[w_j(b)]) / (1 + beta*var2[w_j(k)]),

    where ``w_j`` is the permutation just before the step and ``n(J)`` counts the
    level-``1`` steps.
    """
    chain = epsilon_chain(k, n)
    terms = {}

    def walk(start, w, character, sign):
        terms[w] = terms.get(w, S.Zero) + sign * character / (S.One + beta * var2[w[k - 1]])
        for index in range(start, len(chain)):
            a, b, level = chain[index]
            stepped = w.swap(a - 1, b - 1)
            if stepped.inv != w.inv + 1:
                continue
            if level:
                # r = s_{alpha_kb, 1} translates the weight by w(alpha_kb).
                walk(index + 1, stepped, character * (S.One + beta * var2[w[b - 1]]) / (S.One + beta * var2[w[k - 1]]), sign * beta)
            else:
                walk(index + 1, stepped, character, -sign * beta)

    walk(0, u, S.One, S.One)
    return terms


def single_variable_groth(coeff_dict, varnum, var2=None, beta=None, n=None):
    r"""Multiply ``sum_u coeff_u G_u(x, var2)`` by the single variable ``x_varnum``.

    Returns ``{w: coeff_w}``.  This is ``_one_plus_beta_x_terms`` with the
    identity subtracted off and ``beta`` divided out; the diagonal coefficient
    collapses to ``(-) var2[u(varnum)] = -y/(1 + beta*y)``, the formal inverse of
    ``var2[u(varnum)]``, which is the localization of ``x_varnum`` at ``u``.
    """
    if beta is None:
        beta = _default_beta
    var2 = _genset(var2)
    k = varnum

    ret = {}
    for u, val in coeff_dict.items():
        u = Permutation(u)
        rank = max(n or 0, len(u), k + 1) + 1
        for w, coeff in _one_plus_beta_x_terms(u, k, var2, beta, rank).items():
            if w == u:
                y = var2[u[k - 1]]
                ret[w] = ret.get(w, S.Zero) + val * (-y) / (S.One + beta * y)
            else:
                ret[w] = ret.get(w, S.Zero) + val * _divide_by_beta(coeff, beta)
    return ret


def mult_poly_groth_double(coeff_dict, poly, var_x=None, var_y=None, beta=None, n=None):
    """Multiply ``sum_u coeff_u G_u(x, var_y)`` by an arbitrary polynomial in ``var_x``.

    Mirrors ``mult_poly_double``; the leaves of the ``Add``/``Mul``/``Pow`` recursion
    are handled by ``single_variable_groth``.
    """
    var_x = _genset(var_x)
    var_y = _genset(var_y)
    index = var_x.index(poly)
    if index != -1:
        return single_variable_groth(coeff_dict, index, var_y, beta, n)
    if isinstance(poly, Mul):
        ret = coeff_dict
        for arg in poly.args:
            ret = mult_poly_groth_double(ret, arg, var_x, var_y, beta, n)
        return ret
    if isinstance(poly, Pow):
        base, exponent = poly.args
        ret = coeff_dict
        for _ in range(int(exponent)):
            ret = mult_poly_groth_double(ret, base, var_x, var_y, beta, n)
        return ret
    if isinstance(poly, Add):
        ret = {}
        for arg in poly.args:
            ret = add_perm_dict(ret, mult_poly_groth_double(coeff_dict, arg, var_x, var_y, beta, n))
        return ret
    return {perm: poly * coeff for perm, coeff in coeff_dict.items()}


def _chain_sums(u, k, n, beta):
    """``{w: sum_J beta**(|J| - 1)}`` over nonempty ``J`` with ``u r_J = w``.

    ``J`` runs over subsets of ``monk_chain(k, n)`` whose reflections, applied in
    chain order, form a saturated increasing chain in Bruhat order from ``u``.
    """
    chain = monk_chain(k, n)
    sums = {}

    def walk(start, w, size):
        if size:
            sums[w] = sums.get(w, S.Zero) + beta ** (size - 1)
        for index in range(start, len(chain)):
            a, b = chain[index]
            stepped = w.swap(a - 1, b - 1)
            if stepped.inv == w.inv + 1:
                walk(index + 1, stepped, size + 1)

    walk(0, u, 0)
    return sums


def _genset(var):
    if var is None or isinstance(var, GeneratingSet_base):
        return var
    return CustomGeneratingSet(var)


def _divide_by_beta(expr, beta):
    from sympy import cancel

    return sympify(cancel(sympify_sympy(expr) / sympify_sympy(beta)))


def grothmult_double(perm_dict, v, var2=None, var3=None, beta=None, n=None):
    r"""Multiply double Grothendieck polynomials, mirroring ``schubmult_double``.

    Computes the expansion of ``sum_u coeff_u G_u(x, var2) * G_v(x, var3)`` in
    the basis ``{G_w(x, var2)}`` and returns it as ``{w: coeff_w}``.

    Only ``v.inv <= 1`` is implemented: ``v = 1`` is the identity and ``v = s_k``
    is the degree-one (plus higher ``beta``-corrected) Grothendieck polynomial of
    Corollary 8.2.

    ``n`` overrides the ambient rank used for the ``(-omega_k)``-chain; the
    default is large enough for the product to stabilize.
    """
    if beta is None:
        beta = _default_beta
    var2 = _genset(var2)
    var3 = _genset(var3)

    v = Permutation(v)
    perm_dict = {Permutation(key): value for key, value in perm_dict.items()}
    if v.inv == 0:
        return perm_dict
    if v.inv > 1:
        if v.max_descent == 1:
            from schubmult.combinatorics.permutation import uncode
            numtimes = v.trimcode[0]
            for index in range(numtimes):
                perm_dict = grothmult_double(perm_dict, uncode([1]), var2, var3[index:], beta, n)
            return perm_dict
        raise NotImplementedError(f"grothmult_double is only implemented for v.inv <= 1 or v.max_descent == 1, got {list(v)} with {v.inv} inversions")
    # v = s_k, so its Lehmer code is [0, ..., 0, 1] with the 1 in position k.
    k = len(v.trimcode)

    ret = {}
    for u, val in perm_dict.items():
        rank = max(n or 0, len(u), k + 1) + k

        numerator = S.One
        denominator = S.One
        for i in range(1, k + 1):
            numerator *= S.One + beta * var3[i]
            denominator *= S.One + beta * var2[u[i - 1]]
        theta = numerator / denominator

        # J = {} contributes (Theta_u - 1)/beta, which is regular at beta = 0.
        diagonal = _divide_by_beta(numerator - denominator, beta) / denominator
        ret[u] = ret.get(u, S.Zero) + val * diagonal

        for w, coeff in _chain_sums(u, k, rank, beta).items():
            ret[w] = ret.get(w, S.Zero) + val * theta * coeff

    return ret
