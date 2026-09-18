<a id="schubmult.rings.combinatorial.rc_graph_ring"></a>

# schubmult.rings.combinatorial.rc\_graph\_ring

`RCGraphRing`: the ring whose basis elements are `RCGraph`s.

Two products live here. ``*`` is the dual (stacking) product from `RCGraph.product`,
defined for every pair and compatible with the Schubert-polynomial coproduct
(``vertical_coproduct``). ``%`` (``rc_product``) is the polynomial product,
currently implemented only when one factor is Grassmannian with a sufficiently
large descent (or in two-row cases via squash decomposition). Crystal operators
and most `RCGraph` methods extend linearly to ring elements (see ``broadcast``).
``schub(perm, n)`` is the sum of all RC graphs of ``perm`` with ``n`` rows, i.e. the
Schubert polynomial as a ring element.

`GrassRCGraphRing` restricts to Grassmannian RC graphs (single descent at the last row).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement"></a>

## RCGraphRingElement Objects

```python
class RCGraphRingElement(CrystalGraphRingElement, SchubertMonomialRingElement)
```

RCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for %; the approach taken is to define the ambiguous term in the Leibniz formula
rather than compute the polynomial product directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.as_terms"></a>

#### as\_terms

```python
def as_terms()
```

Terms ``coeff * rc`` in dict order; the empty RC graph is kept as a symbol rather than collapsing to 1.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by RC graph (sympy printing hook).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.vex"></a>

#### vex

```python
@property
def vex()
```

Linear extension of `RCGraph.vex`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.grass"></a>

#### grass

```python
@property
def grass()
```

Linear extension of `RCGraph.grass`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.squash_product"></a>

#### squash\_product

```python
def squash_product(other_rc)
```

Linear extension of RCGraph.squash_product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.vertical_coproduct"></a>

#### vertical\_coproduct

```python
def vertical_coproduct()
```

Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
and induces the mul product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Coproduct via `RCGraphRing.coproduct_on_basis` (currently pattern-restricted; see there).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an RCGraphRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.almosteq"></a>

#### almosteq

```python
def almosteq(other)
```

Whether ``self - other`` has all-zero coefficients.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an RCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an RCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

Linear extension of RCGraph.crystal_reflection:
For each basis RCGraph, apply its crystal_reflection(index) and collect results.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.weight_reflection"></a>

#### weight\_reflection

```python
def weight_reflection(index)
```

Linear extension of `RCGraph.weight_reflection`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.shiftup"></a>

#### shiftup

```python
def shiftup(k)
```

Linear extension of `RCGraph.shiftup`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.prepend"></a>

#### prepend

```python
def prepend(k)
```

Linear extension of `RCGraph.prepend`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
def zero_out_last_row()
```

Linear extension of `RCGraph.zero_out_last_row`, dropping terms whose last row is nonempty.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.resize"></a>

#### resize

```python
def resize(n)
```

Linear extension of `RCGraph.resize`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.clip"></a>

#### clip

```python
def clip(n)
```

Keep the first ``n`` rows of each RC graph (left factor of `RCGraph.vertical_cut`).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.transpose"></a>

#### transpose

```python
def transpose(length)
```

Linear extension of `RCGraph.transpose`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.project"></a>

#### project

```python
def project()
```

Round-trip through the free algebra: ``from_free_algebra_element(to_free_algebra_element())``
(projects onto the sum-of-all-RC-graphs representatives).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.trim_operator"></a>

#### trim\_operator

```python
def trim_operator(i)
```

Linear extension of `RCGraphRing.trim_operator`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.double_elem_sym_squash"></a>

#### double\_elem\_sym\_squash

```python
def double_elem_sym_squash(weight, yvars, zvars)
```

Linear extension of `RCGraph.double_elem_sym_squash`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.full_double_elem_sym_squash"></a>

#### full\_double\_elem\_sym\_squash

```python
def full_double_elem_sym_squash(p, yvars, zvars)
```

Linear extension of `RCGraph.full_double_elem_sym_squash`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.grass_coaction"></a>

#### grass\_coaction

```python
def grass_coaction()
```

Linear extension of `RCGraphRing.grass_coaction`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.broadcast"></a>

#### broadcast

```python
@property
def broadcast()
```

Proxy that lifts any `RCGraph` method linearly: ``elem.broadcast.method(*args)`` applies
``rc.method(*args)`` to each basis RC graph and sums the results with coefficients.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing"></a>

## RCGraphRing Objects

```python
class RCGraphRing(SchubertMonomialRing, CrystalGraphRing)
```

The ring of `RCGraph`s; see the module docstring. Instances are distinct (hashed by an id counter).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

Wrap an `RCGraph` as a basis element (or re-parent an existing element).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.new"></a>

#### new

```python
def new(x)
```

The basis element for RC graph ``x`` with coefficient 1.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.monomial"></a>

#### monomial

```python
def monomial(*tup)
```

The ring element for the monomial ``x^tup``: product of one-row RC graphs, split recursively.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.elem_sym"></a>

#### elem\_sym

```python
def elem_sym(descent, weight)
```

Sum of all RC graphs for the elementary-symmetric permutation with the given ``weight`` and descent.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.from_free_algebra_element"></a>

#### from\_free\_algebra\_element

```python
def from_free_algebra_element(elem)
```

Convert a free-algebra element (via the word basis) into a sum of ``monomial`` elements.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the RCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.double_schub"></a>

#### double\_schub

```python
def double_schub(perm, coeff_genset, n=None)
```

Return the RCGraphRing element corresponding to the double Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used), with `coeff_genset`
as the set of variables for the double part.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.trim_operator"></a>

#### trim\_operator

```python
def trim_operator(i, rc)
```

Divided difference at ``i`` followed by pulling out row ``i`` (only for terms where that row emptied).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.grass_coaction"></a>

#### grass\_coaction

```python
def grass_coaction(elem: RCGraph)
```

Coaction of the Grassmannian ring: squash-decompose ``elem`` as ``base * grass``, coproduct
``grass`` in `GrassRCGraphRing`, and reattach ``base`` to the left factors.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(rc)
```

Coproduct of a single RC graph, lifted from the free-algebra Schubert coproduct via
`BoundedRCFactorAlgebra`; only implemented for permutations avoiding ``4132``, ``1432``, ``3142``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.old_coproduct_on_basis"></a>

#### old\_coproduct\_on\_basis

```python
def old_coproduct_on_basis(elem)
```

Earlier recursive coproduct (peel the last row, recurse, correct); superseded by ``coproduct_on_basis``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.rc_product"></a>

#### rc\_product

```python
def rc_product(elem1, elem2)
```

Polynomial product (``%``), bilinear extension of ``rc_single_product``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.rc_single_product"></a>

#### rc\_single\_product

```python
def rc_single_product(u_rc, v_rc)
```

Polynomial product of two RC graphs of equal length: handled via squash decomposition for two
rows, and via ``squash_product``/``left_squash`` when one factor is Grassmannian with a large
enough descent; raises ``NotImplementedError`` otherwise.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.potential_products"></a>

#### potential\_products

```python
@cache
def potential_products(left, right, length)
```

Candidate RC graphs that could appear in ``left % right``, obtained by multiplying the
transposes with ``*`` and transposing back.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.potential_prodperms"></a>

#### potential\_prodperms

```python
def potential_prodperms(left, right, length)
```

Permutations of the ``potential_products``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.GrassRCGraphRing"></a>

## GrassRCGraphRing Objects

```python
class GrassRCGraphRing(RCGraphRing)
```

`RCGraphRing` restricted to Grassmannian RC graphs (``_is_valid_grass``); non-Grassmannian terms
are filtered out of sums and products, and ``coproduct_on_basis`` is computed recursively by
vertical cuts.

<a id="schubmult.rings.combinatorial.rc_graph_ring.GrassRCGraphRing.mul_pair"></a>

#### mul\_pair

```python
def mul_pair(a, b)
```

Product of two Grassmannian RC graphs: stack ``b`` (shifted) below ``a`` and keep valid Grassmannian results.

<a id="schubmult.rings.combinatorial.rc_graph_ring.GrassRCGraphRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(elem)
```

Coproduct of a Grassmannian RC graph: one-row case is the standard Pieri splitting; otherwise
cut in half vertically, coproduct each half, and keep pairs whose ``%`` product recovers ``elem``.

