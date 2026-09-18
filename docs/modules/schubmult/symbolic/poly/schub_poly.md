<a id="schubmult.symbolic.poly.schub_poly"></a>

# schubmult.symbolic.poly.schub\_poly

Explicit symbolic formulas for (double) Schubert and Grothendieck polynomials.

The main entry points are `schubpoly` (double Schubert polynomial by the ``pull_out_var``
recursion), `schubpoly_from_elems` (Schubert polynomial as a sum of products of elementary
symmetric polynomials along theta-code v-paths, with a pluggable ``elem_func``),
`grothendieck_poly` (via isobaric divided differences), and the divided-difference operators
`div_diff`/`divide_out_diff`. ``_vars`` holds the default generating sets ``x``, ``y``, ``z``,
``q``. Everything here works on raw SymEngine/SymPy expressions; the ring classes in
`schubmult.rings` call these to expand basis elements.

<a id="schubmult.symbolic.poly.schub_poly.sv_posify"></a>

#### sv\_posify

```python
def sv_posify(val, var2)
```

Rewrite ``val`` in the differences ``var2[i+1] - var2[i]`` of consecutive variables (a
positivity-revealing form), by substituting ``var2[i] = var2[1] + r_1 + ... + r_{i-1}``,
simplifying, and mapping the ``r`` variables back.

<a id="schubmult.symbolic.poly.schub_poly.act"></a>

#### act

```python
def act(w, poly, genset)
```

Permute the variables of ``poly``: ``genset[i] -> genset[w(i)]``.

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_func"></a>

#### elem\_sym\_func

```python
def elem_sym_func(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2)
```

The double elementary symmetric factor attached to one step of the ``schubmult_double`` v-path
recursion: ``e_{k - udiff - vdiff}`` in the ``y`` variables fixed by ``u1 -> u2`` and the ``z``
variables selected by `call_zvars` for ``v1 -> v2``.

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_func_q"></a>

#### elem\_sym\_func\_q

```python
def elem_sym_func_q(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2)
```

Quantum-double variant of `elem_sym_func` (all ``k`` positions of ``u1``/``u2`` are compared).

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_poly_q"></a>

#### elem\_sym\_poly\_q

```python
def elem_sym_poly_q(p, k, varl1, varl2, q_var=_vars.q_var)
```

Quantum double elementary symmetric polynomial ``E_p^q(x_1..x_k; y)``: the usual recursion
plus the term ``q_{k-1} E_{p-2}(x_1..x_{k-2})``.

<a id="schubmult.symbolic.poly.schub_poly.complete_sym_poly"></a>

#### complete\_sym\_poly

```python
def complete_sym_poly(p, k, vrs, vrs2)
```

Factorial complete homogeneous symmetric polynomial ``h_p(vrs[0..k-1] | vrs2)``, computed by
splitting the variable set in half.

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_poly"></a>

#### elem\_sym\_poly

```python
def elem_sym_poly(p, k, varl1, varl2, xstart=0, ystart=0)
```

Factorial elementary symmetric polynomial ``e_p(x_1 - y_1, ..., x_k - y_k)`` style sum over
``varl1[xstart:xstart+k]`` and ``varl2[ystart:]``, computed by a divide-and-conquer split of the
variables (the ``y`` offset shifts by the degree taken from the first half).

<a id="schubmult.symbolic.poly.schub_poly.call_zvars"></a>

#### call\_zvars

```python
@cache
def call_zvars(v1, v2, k, i, min_size=10)
```

Indices of the ``z`` variables entering the elementary symmetric factor for the v-path step
``v1 -> v2`` at position ``i`` with ``k`` variables (cached).

<a id="schubmult.symbolic.poly.schub_poly.q_vector"></a>

#### q\_vector

```python
def q_vector(q_exp, q_var=_vars.q_var)
```

Exponent vector of a monomial in the ``q`` variables (``q_1^a q_2^b -> [a, b]``); ``[]`` for 1,
``None`` if ``q_exp`` is not a ``q`` monomial.

<a id="schubmult.symbolic.poly.schub_poly.monom_sym"></a>

#### monom\_sym

```python
def monom_sym(partition, numvars, genset)
```

Monomial symmetric polynomial ``m_partition(genset[1..numvars])``.

<a id="schubmult.symbolic.poly.schub_poly.xreplace_genvars"></a>

#### xreplace\_genvars

```python
def xreplace_genvars(poly, vars1, vars2)
```

Replace the internal placeholder generating sets ``_vars.var_g1``/``var_g2`` with ``vars1``/``vars2``.

<a id="schubmult.symbolic.poly.schub_poly.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(poly, v1, v2)
```

The quotient ``(poly - poly|_{v1 -> v2}) / (v1 - v2)``, computed structurally on the expression
tree (so it is exact and needs no polynomial division). Objects may override via
``_eval_divide_out_diff``.

<a id="schubmult.symbolic.poly.schub_poly.split_up"></a>

#### split\_up

```python
def split_up(poly, v1, v2)
```

Write ``poly = a + (v1 - v2) * b`` with ``a = poly|_{v1 -> v2}``; returns ``(a, (v1 - v2, b))``.

<a id="schubmult.symbolic.poly.schub_poly.perm_act"></a>

#### perm\_act

```python
def perm_act(val, i, var2=None)
```

Swap ``var2[i]`` and ``var2[i+1]`` in ``val`` (the simple transposition ``s_i`` acting on variables).

<a id="schubmult.symbolic.poly.schub_poly.elem_func_func"></a>

#### elem\_func\_func

```python
def elem_func_func(k, i, v1, v2, vdiff, varl1, varl2, elem_func)
```

Single-sided version of `elem_sym_func` with a pluggable ``elem_func(p, k, xvars, zvars)``,
used by `schubpoly_from_elems`.

<a id="schubmult.symbolic.poly.schub_poly.elem_func_func_mul"></a>

#### elem\_func\_func\_mul

```python
def elem_func_func_mul(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2,
                       elem_func)
```

`elem_sym_func` with a pluggable ``elem_func`` in place of `elem_sym_poly`.

<a id="schubmult.symbolic.poly.schub_poly.schubpoly_from_elems"></a>

#### schubpoly\_from\_elems

```python
def schubpoly_from_elems(v, var_x=None, var_y=None, elem_func=None, mumu=None)
```

Schubert polynomial of ``v`` as a sum over v-paths of products of ``elem_func`` factors.

Uses the strict theta code of ``v^{-1}`` (or the code of the dominant ``mumu`` if given) and
the v-path dictionaries of `schubmult.utils.schub_lib.compute_vpathdicts`; each step
contributes ``elem_func(p, k, xvars, zvars)``. With ``elem_func = elem_sym_poly`` this is the
double Schubert polynomial; other choices give the SEM-basis expansion or, as in
`SchubertBasis.transition_word`, an encoding of the factors.

<a id="schubmult.symbolic.poly.schub_poly.schubpoly_classical_from_elems"></a>

#### schubpoly\_classical\_from\_elems

```python
def schubpoly_classical_from_elems(v, var_x=None, var_y=None, elem_func=None)
```

`schubpoly_from_elems` using the ordinary (non-strict) theta code of ``v^{-1}``.

<a id="schubmult.symbolic.poly.schub_poly.schubpoly"></a>

#### schubpoly

```python
def schubpoly(v, var2=None, var3=None, start_var=1)
```

Double Schubert polynomial ``S_v(var2; var3)`` by recursion on the last descent: pull out the
variable ``var2[n]`` (``n`` the last descent) via ``pull_out_var``, multiplying by factors
``(var2[n] - var3[p])``.

<a id="schubmult.symbolic.poly.schub_poly.div_diff"></a>

#### div\_diff

```python
def div_diff(poly, v1, v2)
```

Divided difference ``(poly - s(poly)) / (v1 - v2)`` where ``s`` swaps ``v1`` and ``v2``, computed
structurally on the expression tree. Objects may override via ``_eval_div_diff``.

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly_legacy"></a>

#### grothendieck\_poly\_legacy

```python
@cache
def grothendieck_poly_legacy(perm, x, y, beta, keep_as_schub=False)
```

Double Grothendieck polynomial by descending from ``w0`` (product of ``x (+) y`` factors) with
isobaric divided differences. Superseded by `grothendieck_poly`.

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly"></a>

#### grothendieck\_poly

```python
@cache
def grothendieck_poly(perm, x, y, beta, keep_as_schub=False)
```

Double Grothendieck polynomial ``G_perm(x; y)`` with parameter ``beta``, as an expression or
(``keep_as_schub``) as its double Schubert expansion. See `grothendieck_poly_with_ring`.

<a id="schubmult.symbolic.poly.schub_poly.dom_groth"></a>

#### dom\_groth

```python
@cache
def dom_groth(dom_perm, ring, beta)
```

Double Schubert expansion of the Grothendieck polynomial of a dominant permutation: builds
the product of factorial elementary symmetric factors row by row (from the code of
``dom_perm^{-1}``) with the ``1 + beta y`` twists.

<a id="schubmult.symbolic.poly.schub_poly.isobaric_strip_on_dschub_dict"></a>

#### isobaric\_strip\_on\_dschub\_dict

```python
def isobaric_strip_on_dschub_dict(start, length, perm_dict, coeff_genset,
                                  beta)
```

Apply one isobaric strip to a whole ``{perm: coeff}`` dict, folded.

Coefficients landing on the same permutation merge at every stage instead of
being carried per input basis element, mirroring ``compute_vpathdicts``.

<a id="schubmult.symbolic.poly.schub_poly.isobaric_strip_on_dschub"></a>

#### isobaric\_strip\_on\_dschub

```python
def isobaric_strip_on_dschub(start, length, schub_perm, ring, beta)
```

`isobaric_strip_on_dschub_dict` on a single basis element, returned as a ring element.

<a id="schubmult.symbolic.poly.schub_poly.apply_isobaric_to_schub_dict"></a>

#### apply\_isobaric\_to\_schub\_dict

```python
def apply_isobaric_to_schub_dict(diff_perm, perm_dict, coeff_genset, beta)
```

Fold every strip of ``diff_perm`` over the whole dict, merging between strips.

<a id="schubmult.symbolic.poly.schub_poly.apply_isobaric_to_schub"></a>

#### apply\_isobaric\_to\_schub

```python
@cache
def apply_isobaric_to_schub(diff_perm, schub_perm, ring, beta)
```

`apply_isobaric_to_schub_dict` on a single basis element, returned as a ring element.

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly_with_ring"></a>

#### grothendieck\_poly\_with\_ring

```python
@cache
def grothendieck_poly_with_ring(perm, ring, beta, keep_as_schub=False)
```

Double Grothendieck polynomial via the minimal dominant permutation above ``perm``: start
from `dom_groth` and apply the isobaric divided differences of ``perm^{-1} * dom_perm`` strip
by strip (`apply_isobaric_to_schub_dict`).

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly2"></a>

#### grothendieck\_poly2

```python
@cache
def grothendieck_poly2(perm, x, y, beta, keep_as_schub=False)
```

Variant of `grothendieck_poly_legacy` with ``x - y - beta x y`` factors for ``w0``.

<a id="schubmult.symbolic.poly.schub_poly.to_groth"></a>

#### to\_groth

```python
def to_groth(val, x, y, beta)
```

Expand a polynomial in the double Grothendieck basis ``{perm: coeff}`` by triangular
elimination on monomials: peel off the lowest monomial ``x^c`` (lowest total degree, then lex),
subtract ``coeff * G_{uncode(c)}``, and recurse.

<a id="schubmult.symbolic.poly.schub_poly.to_groth_with_ring"></a>

#### to\_groth\_with\_ring

```python
def to_groth_with_ring(_val, ring, beta)
```

Expand a double Schubert ring element in the double Grothendieck basis.

Triangular elimination by length: for the smallest remaining permutation ``w``, apply the
isobaric divided differences of ``w`` and evaluate at ``x_i = -y_i / (1 + beta y_i)`` (the
point where all nontrivial Grothendieck polynomials vanish) to read off the coefficient of
``G_w``, then subtract ``coeff * G_w`` and repeat. Coefficients are simplified with SymPy.

<a id="schubmult.symbolic.poly.schub_poly.to_groth_with_ring_functional"></a>

#### to\_groth\_with\_ring\_functional

```python
def to_groth_with_ring_functional(_val, ring, beta)
```

`to_groth_with_ring` using the ring element's own ``isobaric_perm`` method.

<a id="schubmult.symbolic.poly.schub_poly.groth_dict_to_poly"></a>

#### groth\_dict\_to\_poly

```python
def groth_dict_to_poly(groth_dict, x, zz, beta)
```

Sum ``coeff * G_perm(x; zz)`` over a ``{perm: coeff}`` dict.

<a id="schubmult.symbolic.poly.schub_poly.schub_elem_to_groth_elem_dict"></a>

#### schub\_elem\_to\_groth\_elem\_dict

```python
@cache
def schub_elem_to_groth_elem_dict(the_perm, beta)
```

Signed count, by ``(inv, max_descent)``, of the permutations ``co_pipe_dream(rc).perm * w0`` over
RC graphs of ``the_perm``, weighted ``(-beta)^(inv difference)``: the Grothendieck-side image of a
Schubert basis element.

<a id="schubmult.symbolic.poly.schub_poly.schub_elem_sym_to_groth_elem_sym_dict"></a>

#### schub\_elem\_sym\_to\_groth\_elem\_sym\_dict

```python
@cache
def schub_elem_sym_to_groth_elem_sym_dict(p, k, beta)
```

`schub_elem_to_groth_elem_dict` for the Grassmannian permutation of ``e_p(x_1..x_k)``, i.e. the
expansion of the elementary symmetric polynomial into Grothendieck-Pieri pieces ``(inv, numvars)``.

<a id="schubmult.symbolic.poly.schub_poly.isobar_it"></a>

#### isobar\_it

```python
def isobar_it(i, genset, elem)
```

K-theoretic isobaric operator ``pi_i`` on a Schubert element: ``partial_i((1 + x_{i+1}) x_i * elem)``
via the nil-Hecke ring (``beta = 1``).

<a id="schubmult.symbolic.poly.schub_poly.lascoux_poly"></a>

#### lascoux\_poly

```python
def lascoux_poly(composition, genset)
```

Lascoux polynomial of a weak composition (``beta = 1``), expanded.

<a id="schubmult.symbolic.poly.schub_poly.groth_elem_as_schub_dict"></a>

#### groth\_elem\_as\_schub\_dict

```python
@cache
def groth_elem_as_schub_dict(perm, beta)
```

Schubert expansion ``{perm': coeff}`` of the Grothendieck polynomial ``G_perm`` (via
``WCGraph.groth_to_schub``).

<a id="schubmult.symbolic.poly.schub_poly.groth_mul_full"></a>

#### groth\_mul\_full

```python
def groth_mul_full(perm_dict, p2, _x, _zz, beta)
```

Multiply a Grothendieck expansion ``perm_dict`` by ``G_p2``: expand ``G_p2`` in Schubert
polynomials and push each through `schub_dict_to_groth_dict`.

<a id="schubmult.symbolic.poly.schub_poly.groth_mul_full_with_ring"></a>

#### groth\_mul\_full\_with\_ring

```python
def groth_mul_full_with_ring(perm_dict, p2, ring, beta)
```

`groth_mul_full` using the ring-aware `schub_dict_to_groth_dict_with_ring`.

<a id="schubmult.symbolic.poly.schub_poly.schub_dict_to_groth_dict"></a>

#### schub\_dict\_to\_groth\_dict

```python
def schub_dict_to_groth_dict(base_groth, schub_dict, beta)
```

Multiply the Grothendieck expansion ``base_groth`` by the Schubert polynomial with expansion
``schub_dict``, returning a Grothendieck expansion.

Writes the Schubert polynomial in the CEM (elementary symmetric) basis, converts each
``e_p(x_1..x_k)`` factor to Grothendieck-Pieri pieces with
`schub_elem_sym_to_groth_elem_sym_dict`, and applies ``groth_pieri_mul`` factor by factor.

<a id="schubmult.symbolic.poly.schub_poly.schub_dict_to_groth_dict_with_ring"></a>

#### schub\_dict\_to\_groth\_dict\_with\_ring

```python
def schub_dict_to_groth_dict_with_ring(base_groth, schub_dict, ring, beta)
```

`schub_dict_to_groth_dict` for a specific ``ring`` (uses ``ring.in_CEM_basis`` and
``ring.is_elem_mul_type`` to recognize elementary symmetric factors).

