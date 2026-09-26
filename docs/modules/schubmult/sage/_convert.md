<a id="schubmult.sage._convert"></a>

# schubmult.sage.\_convert

Conversions between schubmult's SymEngine expressions and Sage ring elements.

Index conventions differ by one: schubmult's generating sets are 1-indexed (``x_1`` is the first
variable, ``x_0`` unused) while Sage's Schubert polynomials expand into ``x0, x1, ...``. Every
conversion here shifts accordingly, so ``y_3`` on the schubmult side is ``y_2`` / ``y2`` on the Sage side.

<a id="schubmult.sage._convert.parse_sage_name"></a>

#### parse\_sage\_name

```python
def parse_sage_name(name)
```

``'x3'`` or ``'y_3'`` -> ``('x', 3)`` / ``('y', 3)`` (0-based Sage index); ``None`` if not indexed.

<a id="schubmult.sage._convert.sage_coefficient_to_symengine"></a>

#### sage\_coefficient\_to\_symengine

```python
def sage_coefficient_to_symengine(c)
```

Base-ring scalar (integer or rational) to a SymEngine number.

<a id="schubmult.sage._convert.sage_polynomial_to_symengine"></a>

#### sage\_polynomial\_to\_symengine

```python
def sage_polynomial_to_symengine(p, gensets, named=None)
```

Sage polynomial (finite or infinite polynomial ring) -> SymEngine expression.

``gensets`` maps a letter to a schubmult ``GeneratingSet``; the Sage variable ``a<i>``/``a_<i>``
becomes ``gensets[a][i + 1]``. Letters not in ``gensets`` are created on the fly. ``named`` maps
unindexed Sage variable names (``'beta'``) to SymEngine symbols.

<a id="schubmult.sage._convert.symengine_to_base_ring"></a>

#### symengine\_to\_base\_ring

```python
def symengine_to_base_ring(exprs, B, named=None)
```

Batch-convert SymEngine expressions to elements of the coefficient ring ``B``.

``B`` may be an ``InfinitePolynomialRing``, the fraction field of one (double Grothendieck
coefficients are rational in ``beta`` and ``y``), or a finite ``PolynomialRing`` (e.g. ``R[beta]``).
``named`` maps unindexed SymEngine symbol names to elements of the finite ring underneath ``B``
(``{'β': beta}``).

A generic tree walk doing every ``+``/``*`` through ``InfinitePolynomial`` arithmetic is Python-level
and ~10x slower than the libsingular ring underneath; here the trees are evaluated directly in the
underlying finite ring (grown first to cover every index that occurs) and wrapped without
conversion. Coefficients are never symbolically expanded on the schubmult side; the normal form is
computed by libsingular. Coefficients of one product share most of their subtrees (the same
``(y_i - z_j)`` factors and partial products recur), so values are memoized on the SymEngine node
across the whole batch.

<a id="schubmult.sage._convert.symengine_to_sage"></a>

#### symengine\_to\_sage

```python
def symengine_to_sage(expr, variable, scalar, named=None)
```

SymEngine expression -> Sage element.

``variable(letter, i)`` returns the Sage element for the schubmult symbol ``letter_i`` (1-based ``i``);
``scalar(n)`` converts a Python ``int``/``Fraction``-like rational to the target ring; ``named`` maps
unindexed symbol names (``'β'``) to target elements.

