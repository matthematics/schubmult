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
def sage_polynomial_to_symengine(p, gensets)
```

Sage polynomial (finite or infinite polynomial ring) -> SymEngine expression.

``gensets`` maps a letter to a schubmult ``GeneratingSet``; the Sage variable ``a<i>``/``a_<i>``
becomes ``gensets[a][i + 1]``. Letters not in ``gensets`` are created on the fly.

<a id="schubmult.sage._convert.symengine_to_sage"></a>

#### symengine\_to\_sage

```python
def symengine_to_sage(expr, variable, scalar)
```

SymEngine expression -> Sage element.

``variable(letter, i)`` returns the Sage element for the schubmult symbol ``letter_i`` (1-based ``i``);
``scalar(n)`` converts a Python ``int``/``Fraction``-like rational to the target ring.

