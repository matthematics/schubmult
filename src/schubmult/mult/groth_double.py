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
from schubmult.symbolic import S, sympify, sympify_sympy
from schubmult.symbolic.poly.variables import CustomGeneratingSet, GeneratingSet_base

__all__ = ["grothmult_double", "monk_chain"]


def monk_chain(k, n):
    """Reduced ``(-omega_k)``-chain of reflections in ``A_{n-1}`` (Cor. 15.4).

    Returned as the ordered tuple of transpositions ``(i, j)``, ``i <= k < j``.
    """
    return tuple((i, j) for i in range(1, k + 1) for j in range(n, k, -1))


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
