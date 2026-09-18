<a id="schubmult.rings.thompson_algebra"></a>

# schubmult.rings.thompson\_algebra

`ThompsonAlgebra`: a noncommutative algebra on words in generators ``T_i`` and ``R_i``.

A monomial is a tuple of nonzero integers: ``i > 0`` stands for ``T_i`` and ``i < 0`` for
``R_{-i}``. Products are rewritten to a normal form by the rules in `ThompsonAlgebra._commute_pair`:
``T_i T_j = T_j T_{i+1}`` for ``i > j`` (the Thompson monoid relation), and ``T_i R_j`` moves
``R`` to the left with the index shifts (and one two-term case) given there.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebraElement"></a>

## ThompsonAlgebraElement Objects

```python
class ThompsonAlgebraElement(BaseRingElement)
```

Element of `ThompsonAlgebra`: a dict from normal-form words to coefficients.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra"></a>

## ThompsonAlgebra Objects

```python
class ThompsonAlgebra(BaseRing)
```

Algebra on words in ``T_i`` (positive index) and ``R_i`` (negative index). See the module docstring.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra.printing_term"></a>

#### printing\_term

```python
@cache
def printing_term(monomial)
```

Noncommutative product of the ``T_i``/``R_i`` symbols for the word.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra.new"></a>

#### new

```python
def new(x)
```

Build an element from a word (normalized via `_mul_monomials`), a number, or an existing element.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Scalar multiplication, or the bilinear extension of `_mul_monomials`.

