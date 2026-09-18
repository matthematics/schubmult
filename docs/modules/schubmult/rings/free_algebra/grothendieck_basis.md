<a id="schubmult.rings.free_algebra.grothendieck_basis"></a>

# schubmult.rings.free\_algebra.grothendieck\_basis

`GrothendieckBasis`: the free-algebra basis dual to Grothendieck polynomials.

Keys are ``(perm, numvars)`` as in `SchubertBasis`; the key is dual to ``G_perm`` in ``numvars``
variables (``GrothendieckPolyBasis`` on the polynomial side). Grothendieck polynomials are
taken at ``beta = 1``, which is without loss of generality: ``beta`` is recovered by grading
(a term of ``G_w`` in degree ``inv(w) + d`` carries ``beta^d``). The change of
basis to `SchubertBasis` is the transpose of the Grothendieck-to-Schubert expansion and is
computed combinatorially: enumerate unreduced BPDs of ``perm * w0``, take the co-BPD, and keep
the reduced ones, with sign ``(-1)^(inv(perm) - inv(result))``. Products and all other
transitions route through `SchubertBasis`. ``AGx`` is the standard instance.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis"></a>

## GrothendieckBasis Objects

```python
class GrothendieckBasis(FreeAlgebraBasis)
```

Free-algebra basis dual to Grothendieck polynomials (``beta = 1``); keys are ``(Permutation, numvars)``.

See the module docstring. Products and transitions go through `SchubertBasis`.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Whether ``x`` is ``(perm,)`` or ``(perm, numvars)``.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize to ``(Permutation, numvars)``; ``numvars`` defaults to the last descent.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
@cache
def transition_schubert(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in `SchubertBasis`.

For each unreduced BPD of ``perm * w0`` whose co-BPD is reduced with permutation ``u``
fitting in ``numvars`` variables, contributes ``(-1)^(inv(perm) - inv(u))`` to ``(u, numvars)``.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Key -> ``{key: coeff}`` function into ``other_basis``: identity on Grothendieck subclasses,
`transition_schubert` for `SchubertBasis`, and Schubert-then-onward for everything else.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Display as ``AGx(perm, numvars)``.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply by expanding both keys in `SchubertBasis`, using its separated-descents product,
and converting the result back.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

``GrothendieckPolyBasis``: Grothendieck polynomials are the dual basis.

