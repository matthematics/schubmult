<a id="schubmult.sage._common"></a>

# schubmult.sage.\_common

Shared machinery for Sage parents backed by schubmult rings.

A backed ring is a :class:`~sage.combinat.free_module.CombinatorialFreeModule` indexed by
permutations whose base ring is an infinite polynomial ring over the scalars in the coefficient
alphabets (``y``, ``z``, ``q``, ...). Arithmetic is delegated to a schubmult ring object
(:meth:`SchubmultBackedRing._schub_ring`) and coefficients are converted at the boundary.

Index conventions: Sage variables are 0-indexed (``x0``, ``y_0``, ``q_0``), schubmult's are 1-indexed
(``x_1``, ``y_1``, ``q_1``); see :mod:`schubmult.sage._convert`.

<a id="schubmult.sage._common.coefficient_into"></a>

#### coefficient\_into

```python
def coefficient_into(c, T)
```

Move a base-ring coefficient (infinite polynomial in ``a_<i>``) into the finite ring ``T`` with variables ``a<i>``.

<a id="schubmult.sage._common.SchubmultBackedElement"></a>

## SchubmultBackedElement Objects

```python
class SchubmultBackedElement(CombinatorialFreeModule.Element)
```

<a id="schubmult.sage._common.SchubmultBackedElement.expand"></a>

#### expand

```python
def expand()
```

Expand into a polynomial in ``x0, x1, ...`` and the coefficient variables ``y0, q0, ...``.

There are `n` variables ``x``, `n` the size of the largest permutation involved (as for
:meth:`sage.combinat.schubert_polynomial.SchubertPolynomial_class.expand`); coefficient
letters get exactly the indices that occur.

<a id="schubmult.sage._common.SchubmultBackedRing"></a>

## SchubmultBackedRing Objects

```python
class SchubmultBackedRing(CombinatorialFreeModule)
```

Base class: subclasses set ``_alphabet`` (basis alphabet letter or ``None``), ``_alphabets`` (all
coefficient letters, sorted), and implement ``_schub_ring(alphabet)`` and ``_check_basis_perm``.

