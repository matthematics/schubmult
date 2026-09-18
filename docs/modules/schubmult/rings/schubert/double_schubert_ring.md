<a id="schubmult.rings.schubert.double_schubert_ring"></a>

# schubmult.rings.schubert.double\_schubert\_ring

Double Schubert polynomial ring: the ``DSx`` interface.

`DoubleSchubertRing` represents ``Z[y][x]`` in the basis of double Schubert
polynomials ``S_w(x; y)``, dispatching products to `schubmult.mult.double`.
It is also the workhorse behind the single ring (`schubert_ring.SingleSchubertRing`
is a `DoubleSchubertRing` with an all-zero coefficient alphabet). Beyond ring
arithmetic, `DoubleSchubertElement` supports divided differences, isobaric
divided differences, variable substitution/evaluation, coproducts, and
expansion into elementary-symmetric ("CEM"/"SEM") bases.

Variants: `ElemDoubleSchubertRing` keeps coefficients as unevaluated factorial
elementary symmetric functions; `DoubleSchubertRingDown` uses the descent-side
("down") kernels.

<a id="schubmult.rings.schubert.double_schubert_ring.is_fact_elem_sym"></a>

#### is\_fact\_elem\_sym

```python
def is_fact_elem_sym(obj)
```

Whether ``obj`` is an (unevaluated) factorial elementary symmetric function.

<a id="schubmult.rings.schubert.double_schubert_ring.is_fact_complete_sym"></a>

#### is\_fact\_complete\_sym

```python
def is_fact_complete_sym(obj)
```

Whether ``obj`` is an (unevaluated) factorial complete homogeneous symmetric function.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement"></a>

## DoubleSchubertElement Objects

```python
class DoubleSchubertElement(BaseSchubertElement)
```

An element of a `DoubleSchubertRing`: ``{Permutation: coefficient}`` in the
double Schubert basis ``S_w(x; y)``, with sympy coefficients in ``y``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.to_genset_dict"></a>

#### to\_genset\_dict

```python
def to_genset_dict(trim=False)
```

Expand to a polynomial and return ``{exponent_tuple: coeff}`` over the ``x`` variables;
``trim=True`` merges keys that differ only by trailing zeros.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.divdiff"></a>

#### divdiff

```python
def divdiff(i)
```

Divided difference ``partial_i``: ``S_w -> S_{w s_i}`` when ``i`` is a descent of ``w``, else 0.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.simpleref"></a>

#### simpleref

```python
def simpleref(i)
```

Action of the simple reflection ``s_i`` on the ``x`` variables: ``f + (x_{i+1} - x_i) partial_i f``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.coeff_isobaric"></a>

#### coeff\_isobaric

```python
def coeff_isobaric(i, beta)
```

Isobaric divided difference acting on the ``y`` (coefficient) alphabet, transported through
the basis via the antipode-style inversion ``S_w -> (-1)^{l(w)} S_{w^{-1}}``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.isobaric"></a>

#### isobaric

```python
def isobaric(i, beta)
```

Beta-deformed isobaric divided difference ``pi_i = partial_i + beta (x_i partial_i - 1)``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply ``partial_w`` for ``w = perm``, peeling simple reflections from the last descent.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.isobaric_perm"></a>

#### isobaric\_perm

```python
def isobaric_perm(perm, beta)
```

Apply the beta-isobaric ``pi_w`` for ``w = perm``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.isobaric_plus_beta"></a>

#### isobaric\_plus\_beta

```python
def isobaric_plus_beta(i, beta)
```

The variant ``partial_i + beta x_i partial_i`` (isobaric without the ``-beta`` identity term).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.act"></a>

#### act

```python
def act(perm)
```

Permute the ``x`` variables by ``perm``, as a composition of ``simpleref``s.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.max_index"></a>

#### max\_index

```python
def max_index()
```

The largest ``x`` index (1-indexed) any basis permutation actually depends on.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.eval"></a>

#### eval

```python
def eval(x)
```

Substitute ``{generator: value}`` pairs one at a time (via ``pull_out_gen``); returns a
scalar if the result collapses to the identity basis element.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.subs"></a>

#### subs

```python
def subs(old, new)
```

Substitute ``old -> new`` where ``old`` is an ``x`` variable (moved to the last position and
pulled out via ``pull_out_var``), a ``y`` variable (transported through the basis), or a plain
coefficient symbol.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.free_symbols"></a>

#### free\_symbols

```python
@property
def free_symbols()
```

Coefficient symbols plus the ``x``/``y`` variables the basis permutations actually depend on.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.pull_out_gen"></a>

#### pull\_out\_gen

```python
def pull_out_gen(gen)
```

Factor out all dependence on one generator ``gen`` (an ``x`` or ``y`` variable), returning an
element over a `MaskedGeneratingSet` ring with ``gen`` removed and explicit ``(gen - y_j)``
(or factorial-elementary-symmetric) prefactors.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.in_CEM_basis"></a>

#### in\_CEM\_basis

```python
def in_CEM_basis()
```

Expand in the complete-elementary-monomial (CEM) basis using the ring's symbolic elementary function.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.cem_rep"></a>

#### cem\_rep

```python
def cem_rep(elem_func, mumu=None)
```

CEM expansion with a custom ``elem_func``; ``mumu`` selects a dominant permutation to expand
against (defaults to the classical route).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.coproduct"></a>

#### coproduct

```python
def coproduct(*indices,
              alt_coeff_genset=None,
              on_coeff_gens=False,
              gname1=None,
              gname2=None)
```

Coproduct splitting the ``x`` variables (or ``y`` if ``on_coeff_gens``) at the given 1-indexed
``indices``: returns an element of the `TensorRing` of two `DoubleSchubertRing`s over the
complementary `MaskedGeneratingSet`s, labeled ``gname1``/``gname2``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.max_gens"></a>

#### max\_gens

```python
@cached_property
def max_gens()
```

Largest 0-indexed descent over all basis permutations.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.positive_elem_sym_rep"></a>

#### positive\_elem\_sym\_rep

```python
def positive_elem_sym_rep()
```

Manifestly positive expansion in factorial elementary symmetric functions (forward ``pull_out_var``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.positive_elem_sym_rep_backward"></a>

#### positive\_elem\_sym\_rep\_backward

```python
def positive_elem_sym_rep_backward()
```

Like ``positive_elem_sym_rep`` but peeling from the last descent backward.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.antipode"></a>

#### antipode

```python
def antipode()
```

The antipode: swap the two alphabets and invert each basis permutation (see `DoubleSchubertRing.antipode`).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing"></a>

## DoubleSchubertRing Objects

```python
class DoubleSchubertRing(BaseSchubertRing)
```

The ring of double Schubert polynomials ``S_w(x; y)`` over ``genset`` (``x``) and
``coeff_genset`` (``y``). Call the ring with a permutation, Lehmer code, or polynomial
expression to construct an element; the module-level ``DSx`` is the standard instance.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.coeff_ring"></a>

#### coeff\_ring

```python
@cached_property
def coeff_ring()
```

The single Schubert ring over the coefficient alphabet ``y`` (used by ``coeff_isobaric``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.antipode_ring"></a>

#### antipode\_ring

```python
@cached_property
def antipode_ring()
```

The same ring with the two alphabets swapped.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.antipode"></a>

#### antipode

```python
def antipode(elem)
```

Map ``sum c_w S_w(x; y)`` to ``sum c_w S_{w^{-1}}(y; x)`` in the swapped-alphabet ring.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Right-multiply by a scalar (coefficient-domain element) or, failing that, by an expression.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.positive_elem_sym_rep"></a>

#### positive\_elem\_sym\_rep

```python
def positive_elem_sym_rep(perm, index=1)
```

Manifestly positive expansion of ``S_perm`` in factorial elementary symmetric functions, peeling
the first variable of ``~perm`` at each step (``pull_out_var(1, ...)``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.positive_elem_sym_rep_backward"></a>

#### positive\_elem\_sym\_rep\_backward

```python
def positive_elem_sym_rep_backward(perm)
```

Like ``positive_elem_sym_rep`` but peeling from the last descent of ``~perm`` backward.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="")
```

The ``DSchubPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
@property
def elem_sym()
```

`FactorialElemSym`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(other)
```

Whether ``other`` is a factorial elementary symmetric function (eligible for ``elem_mul``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Multiply by a factorial elementary symmetric function in ``x`` variables via the positional
Pieri rule (``elem_sym_positional_perms``), expanding the leftover factor with ``expand_func``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

`FactorialElemSym` (kept unevaluated for symbolic expansions).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.schubert_schur_elem_func"></a>

#### schubert\_schur\_elem\_func

```python
def schubert_schur_elem_func(numvars)
```

Elementary-symmetric substitute for the Schubert-tensor-Schur expansion: ``e_p(x_1..x_k)`` maps
to a Schubert basis element on the left factor when ``k >= numvars`` and on the right otherwise.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_schubert_schur_basis"></a>

#### in\_schubert\_schur\_basis

```python
def in_schubert_schur_basis(perm, numvars)
```

Expand ``S_perm`` in the Schubert-tensor-Schur basis, treating the last ``numvars`` variables
as the symmetric (Schur) part.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_descending_schur_basis"></a>

#### in\_descending\_schur\_basis

```python
def in_descending_schur_basis(perm, numvars)
```

Iterate ``in_schubert_schur_basis`` down through ``numvars, numvars-1, ..., 1``, producing a
nested tensor of Schur-like factors.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Substitution dict ``{e_p_k: elem_sym_poly(p, k, x)}`` for all ``1 <= p <= k <= kk``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.flip"></a>

#### flip

```python
@staticmethod
def flip(elem)
```

Re-express a factorial elementary symmetric function with its two alphabets swapped, via the
corresponding Grassmannian Schubert polynomial's CEM expansion.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Expand each basis element via ``quantum_schubpoly`` (a quantum double Schubert element).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Identity (this ring is already classical).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.quantum_schubpoly"></a>

#### quantum\_schubpoly

```python
@cache
def quantum_schubpoly(perm)
```

The classical ``S_perm`` expressed in the quantum double Schubert basis (via ``quantum_elem_func``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants of ``S_u(x; y) * S_v(x; z)`` (``z`` = ``basis2.coeff_genset``), via ``schubmult_double``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Like ``cached_product`` but with manifestly positive coefficients (generic alphabets, then substituted).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.double.schubmult_double`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.single.schubmult_py`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.single.mult_poly_py`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.double.mult_poly_double`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.quantum_elem_func"></a>

#### quantum\_elem\_func

```python
@property
def quantum_elem_func()
```

Elementary symmetric function valued in the quantum double Schubert ring, computed by a
divide-and-conquer recursion on the variable set (used by ``quantum_schubpoly``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.monomial_schub"></a>

#### monomial\_schub

```python
def monomial_schub(monom)
```

The monomial ``x^monom`` expressed in the Schubert basis (trailing zeros in ``monom`` ignored).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit polynomial ``S_k(x; y)`` (cached).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.complete_mul"></a>

#### complete\_mul

```python
def complete_mul(elem, x)
```

Multiply by a factorial complete homogeneous symmetric function in ``x`` variables via
``complete_sym_positional_perms`` (the dual Pieri rule).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.handle_sympoly"></a>

#### handle\_sympoly

```python
def handle_sympoly(other)
```

How a symmetric-function coefficient is stored: evaluated to a polynomial here.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.single_variable"></a>

#### single\_variable

```python
def single_variable(elem, varnum)
```

Multiply by the single variable ``x_varnum`` (equivariant Monk rule).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial expression in ``x``/``y`` into the Schubert basis.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply ``elem`` by an arbitrary expression ``x``: single variables use the Monk rule,
(factorial) elementary/complete symmetric functions use their Pieri rules (splitting out
variables from the wrong alphabet as needed), and ``Add``/``Mul``/``Pow`` recurse; anything
else is treated as a coefficient.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an existing element of this ring, or an expression.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown"></a>

## DoubleSchubertRingDown Objects

```python
class DoubleSchubertRingDown(DoubleSchubertRing)
```

`DoubleSchubertRing` using the descent-side ("down") multiplication kernels
(``schubmult_double_down``/``schubmult_py_down``); basis symbols print with an ``op`` prefix.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.double.schubmult_double_down`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.single.schubmult_py_down`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Down-kernel structure constants over generic alphabets, substituted back to the ring's alphabets.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Positive variant of ``cached_product`` for the down kernel.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="op")
```

The ``DSchubPoly`` display symbol, prefixed with ``op`` by default.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing"></a>

## ElemDoubleSchubertRing Objects

```python
class ElemDoubleSchubertRing(DoubleSchubertRing)
```

`DoubleSchubertRing` whose coefficients are kept as unevaluated `FactorialElemSym`
functions instead of being expanded to polynomials; products use the ``*_from_elems`` kernels.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.replacematch"></a>

#### replacematch

```python
@property
def replacematch()
```

A ``(a, b) -> expression`` rewriter turning differences ``a - b`` into `FactorialElemSym(1, 1, ...)`
forms, respecting which alphabet each symbol belongs to.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.elem_func"></a>

#### elem\_func

```python
@property
def elem_func()
```

`FactorialElemSym`.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.handle_sympoly"></a>

#### handle\_sympoly

```python
def handle_sympoly(other)
```

Keep symmetric-function coefficients unevaluated.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Positional Pieri rule for a factorial elementary symmetric function, keeping the leftover
factor as an unevaluated coefficient.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.complete_mul"></a>

#### complete\_mul

```python
def complete_mul(elem, x)
```

Dual Pieri rule for a factorial complete symmetric function, keeping the leftover factor unevaluated.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants via ``schubmult_double_from_elems`` with `FactorialElemSym` coefficients.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Structure constants via the positive ``schubmult_double_alt_from_elems`` route.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an element of this ring, or an expression.

<a id="schubmult.rings.schubert.double_schubert_ring.DSx"></a>

#### DSx

```python
def DSx(x, genset=GeneratingSet("y"), elem_sym=False, down=False)
```

Construct a double Schubert polynomial element in ``x`` with coefficient alphabet ``genset``.

``DSx([3, 1, 2])`` is ``S_{312}(x; y)``. Pass ``genset="z"`` (or a `GeneratingSet`) for a
different coefficient alphabet; ``elem_sym=True`` uses `ElemDoubleSchubertRing`, ``down=True``
uses `DoubleSchubertRingDown`.

