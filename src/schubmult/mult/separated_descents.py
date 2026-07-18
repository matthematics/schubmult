r"""Separated-descents product of (double) Grothendieck polynomials.

Implements the *pipe puzzle* formula of Fan--Guo--Xiong,
"Bumpless pipe dreams meet puzzles" (arXiv:2309.00467), Theorem 2.5.

For permutations ``u`` and ``v`` with *separated descents* at a position ``k``
(i.e. every descent of ``u`` is ``<= k`` and every descent of ``v`` is ``>= k``)
the double Grothendieck polynomials satisfy

    G_u(x, y) * G_v(x, t) = sum_w c_{u,v}^w(t, y) * G_w(x, t),

and the structure constants ``c_{u,v}^w(t, y)`` are given by a positive
(in the sense of Theorem 2.5) sum over pipe puzzles.

The beta convention matches the rest of ``schubmult``: the multiplicative
formal group law is ``x (+) y = x + y + beta*x*y`` with formal subtraction

    x (-) y = (x - y) / (1 + beta*y).

Setting ``beta = 0`` recovers the double Schubert (cohomology) structure
constants of Theorem 4.x (the "Schubert pipe puzzle" specialization).
"""

from schubmult.abc import beta as _default_beta
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import S, expand, sympify
from schubmult.utils.perm_utils import add_perm_dict

__all__ = ["grothmult_double", "separated_descents_coeffs"]


def _grothendieck_x_degree(perm):
    """Top total degree in the ``x`` variables of the Grothendieck polynomial ``G_perm``.

    Computed with the secondary variables set to zero and ``beta`` specialized
    to ``1`` (so it does not contribute to the degree count).
    """
    import sympy

    from schubmult.abc import x
    from schubmult.symbolic.poly.schub_poly import grothendieck_poly
    from schubmult.symbolic.poly.variables import ZeroGeneratingSet

    if perm.inv == 0:
        return 0
    poly = expand(grothendieck_poly(perm, x, ZeroGeneratingSet(), 1))
    expr = sympy.sympify(str(poly))
    syms = sorted(expr.free_symbols, key=str)
    if not syms:
        return 0
    return sympy.Poly(expr, *syms).total_degree()


def _one_line(perm, n):
    """Return the one-line notation of ``perm`` padded with fixed points to length ``n``."""
    arr = list(perm)
    return arr + list(range(len(arr) + 1, n + 1))


def _inverse_one_line(arr):
    """Inverse of a one-line permutation array (1-indexed values)."""
    n = len(arr)
    inv = [0] * n
    for pos in range(n):
        inv[arr[pos] - 1] = pos + 1
    return inv


def _separating_position(u, v):
    """Return ``k`` witnessing separated descents of ``u`` and ``v``.

    Raises ``ValueError`` if no separating position exists.
    """
    des_u = u.descents(zero_indexed=False)
    des_v = v.descents(zero_indexed=False)
    max_u = max(des_u) if des_u else 0
    min_v = min(des_v) if des_v else None
    if min_v is not None and max_u > min_v:
        raise ValueError(
            f"permutations {list(u)} and {list(v)} do not have separated descents "
            f"(max descent of u is {max_u}, min descent of v is {min_v})",
        )
    return max_u


def _coeffs_on_grid(u, v, var1, var2, beta, k, n):
    """Enumerate pipe puzzles on a fixed ``n`` by ``n`` grid (Theorem 2.5).

    Returns a dict ``{w: c_{u,v}^w}`` (coefficients not yet simplified).
    """
    u_inv = _inverse_one_line(_one_line(u, n))
    v_inv = _inverse_one_line(_one_line(v, n))

    # Right-side boundary labels (rows, top to bottom): kappa_u^i.
    kappa = [u_inv[i] if u_inv[i] <= k else 0 for i in range(n)]
    # Top-side boundary labels (columns, left to right): theta_v^j.
    theta = [v_inv[j] if v_inv[j] > k else 0 for j in range(n)]

    tvar = var1
    yvar = var2

    def ominus(a, b):
        return (a - b) / (S.One + beta * b)

    results = {}

    def record(bottom, weight):
        # ``bottom`` is w^{-1} in one-line notation (column -> label).
        w = ~Permutation(list(bottom))
        results[w] = results.get(w, S.Zero) + weight

    def options(t_in, r_in, at_left, i, j):
        """Yield ``(bottom_out, left_out, factor)`` for the tiles fitting a cell.

        ``i``, ``j`` are 1-indexed row/column (for variable subscripts).
        """
        empty = ominus(tvar[j], yvar[i])
        groth = S.One + beta * empty
        if t_in == 0 and r_in == 0:
            yield (0, 0, empty)  # \O
            return
        if t_in != 0 and r_in == 0:
            yield (t_in, 0, S.One)  # \I
            if not at_left:
                # \J : pipe turns top -> left; nontrivial weight if from right side (label <= k)
                yield (0, t_in, groth if t_in <= k else S.One)
            return
        if t_in == 0 and r_in != 0:
            # \F : pipe turns right -> bottom; nontrivial weight if from top side (label > k)
            yield (r_in, 0, groth if r_in > k else S.One)
            if not at_left:
                yield (0, r_in, S.One)  # \H
            return
        # Both pipes present.
        if at_left:
            return  # neither \X nor \B can send a pipe out the left boundary
        if r_in < t_in:
            # \X cross: horizontal label must be smaller than vertical label
            yield (t_in, r_in, S.One)
        same_side = (t_in <= k) == (r_in <= k)
        if same_side:
            # \B bump, both from same side: NW (top-left, label t_in) must be larger
            if t_in > r_in:
                yield (r_in, t_in, beta)
        elif t_in <= k and r_in > k:
            # \B bump, different sides: NW pipe must enter from the right side
            yield (r_in, t_in, beta * groth)

    def rec(i, j, top, right_label, weight):
        if j < 0:
            if i == n - 1:
                record(top, weight)
                return
            rec(i + 1, n - 1, top, kappa[i + 1], weight)
            return
        t_in = top[j]
        r_in = right_label
        for bottom_out, left_out, factor in options(t_in, r_in, j == 0, i + 1, j + 1):
            saved = top[j]
            top[j] = bottom_out
            rec(i, j - 1, top, left_out, weight * factor)
            top[j] = saved

    rec(0, n - 1, list(theta), kappa[0], S.One)

    return results


def separated_descents_coeffs(u, v, var1, var2, beta=None, grid_size=None):
    r"""Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

    ``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
    ``var2`` are the ``y`` variables (secondary variables of ``u``).

    Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

    In K-theory the product may involve ``G_w`` with ``w`` in a larger
    symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
    ``grid_size`` is not supplied it is chosen large enough to capture every
    such ``w``: the maximum value of any appearing ``w`` is bounded by
    ``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.
    """
    if beta is None:
        beta = _default_beta
    u = Permutation(u)
    v = Permutation(v)
    k = _separating_position(u, v)

    def clean(raw):
        return {w: c for w, c in raw.items() if expand(sympify(c)) != S.Zero}

    if grid_size is None:
        n = max(len(u), len(v), 2)
        grid_size = max(n, (n - 1) + _grothendieck_x_degree(u) + _grothendieck_x_degree(v))

    return clean(_coeffs_on_grid(u, v, var1, var2, beta, k, grid_size))


def grothmult_double(perm_dict, v, var1, var2, beta=None):
    r"""Separated-descents product of double Grothendieck polynomials.

    Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
    ``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

        coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

    with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

    ``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.
    """
    if beta is None:
        beta = _default_beta
    v = Permutation(v)
    ret = {}
    for u, coeff in perm_dict.items():
        single = separated_descents_coeffs(u, v, var1, var2, beta=beta)
        ret = add_perm_dict(ret, {w: c * coeff for w, c in single.items()})
    return {w: c for w, c in ret.items() if expand(sympify(c)) != S.Zero}
