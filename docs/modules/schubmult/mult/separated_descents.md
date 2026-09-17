<a id="schubmult.mult.separated_descents"></a>

# schubmult.mult.separated\_descents

Separated-descents product of (double) Grothendieck polynomials.

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

<a id="schubmult.mult.separated_descents.separated_descents_coeffs"></a>

#### separated\_descents\_coeffs

```python
def separated_descents_coeffs(u, v, var1, var2, beta=None, grid_size=None)
```

Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
``var2`` are the ``y`` variables (secondary variables of ``u``).

Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

In K-theory the product may involve ``G_w`` with ``w`` in a larger
symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
``grid_size`` is not supplied it is chosen large enough to capture every
such ``w``: the maximum value of any appearing ``w`` is bounded by
``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.

<a id="schubmult.mult.separated_descents.separated_descents_coeffs_plus"></a>

#### separated\_descents\_coeffs\_plus

```python
def separated_descents_coeffs_plus(u,
                                   v,
                                   var1,
                                   var2,
                                   beta=None,
                                   grid_size=None,
                                   mangle_genset=False)
```

Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
``var2`` are the ``y`` variables (secondary variables of ``u``).

Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

In K-theory the product may involve ``G_w`` with ``w`` in a larger
symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
``grid_size`` is not supplied it is chosen large enough to capture every
such ``w``: the maximum value of any appearing ``w`` is bounded by
``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.

<a id="schubmult.mult.separated_descents.grothmult_double"></a>

#### grothmult\_double

```python
def grothmult_double(perm_dict, v, var1, var2, beta=None)
```

Separated-descents product of double Grothendieck polynomials.

Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

    coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.

<a id="schubmult.mult.separated_descents.grothmult_double_plus"></a>

#### grothmult\_double\_plus

```python
def grothmult_double_plus(perm_dict,
                          v,
                          var1,
                          var2,
                          beta=None,
                          mangle_genset=False)
```

Separated-descents product of double Grothendieck polynomials.

Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

    coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.

