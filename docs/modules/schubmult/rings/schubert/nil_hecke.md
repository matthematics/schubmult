<a id="schubmult.rings.schubert.nil_hecke"></a>

# schubmult.rings.schubert.nil\_hecke

The nilHecke ring of divided-difference operators acting on Schubert polynomials.

`NilHeckeRing` elements are ``{Permutation: coefficient}`` combinations of the
divided-difference operators ``partial_w`` (printed ``df(w)``), with polynomial
coefficients in the ``x`` variables multiplied on the left. ``partial_w`` acts on
a `DoubleSchubertElement` via `NilHeckeElement.apply`, sending ``S_v -> S_{v w^{-1}}``
when length-additive. Products use the descent-side kernel ``schubmult_double_down``
to commute polynomial coefficients past operators. The module-level ``df`` is the
standard instance in ``x``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement"></a>

## NilHeckeElement Objects

```python
class NilHeckeElement(DomainElement, DefaultPrinting, dict)
```

An element of a `NilHeckeRing`: ``{Permutation: coeff}`` combination of divided-difference operators.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.apply"></a>

#### apply

```python
def apply(other)
```

Act on a `DoubleSchubertElement`: each ``partial_w`` sends ``S_v -> S_{v w^{-1}}`` when
``l(v w^{-1}) = l(v) - l(w)``, else kills it; coefficients multiply the result.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_terms"></a>

#### as\_terms

```python
def as_terms()
```

Terms ``coeff * df(w)`` in dict order (sympy printing hook).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by permutation length then lexicographically (sympy printing hook).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

``{df(w): coeff}`` mapping display symbols to coefficients.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand each coefficient, keeping the operator basis.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_expr"></a>

#### as\_expr

```python
def as_expr()
```

Sum of the ``as_terms()`` as a sympy ``Add``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing"></a>

## NilHeckeRing Objects

```python
class NilHeckeRing(Ring, CompositeDomain)
```

The nilHecke ring in the alphabet ``genset``; see the module docstring. ``df`` is the standard instance.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.to_sympy"></a>

#### to\_sympy

```python
def to_sympy(elem)
```

Convert an element to a sympy expression (``as_expr``).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.isobaric"></a>

#### isobaric

```python
def isobaric(perm, groth=False, *, groth_beta=None)
```

The isobaric divided difference ``pi_perm`` as a nilHecke element: ``pi_i = partial_i x_{i+1}``
(or the Grothendieck version ``partial_i (1 + beta x_{i+1})`` with ``groth=True``), composed
along a reduced word of ``perm``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.g_isobaric"></a>

#### g\_isobaric

```python
def g_isobaric(perm)
```

``isobaric(perm, groth=True)``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.fgp_operator"></a>

#### fgp\_operator

```python
def fgp_operator(k, length, q_var=GeneratingSet("q"))
```

The Fomin-Gelfand-Postnikov quantization of ``x_k`` as a nilHecke element:
``x_k - sum_{i<k} q_i...q_{k-1} partial_{(i k)} + sum_{i>k} q_k...q_{i-1} partial_{(k i)}``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.subs_fgp"></a>

#### subs\_fgp

```python
def subs_fgp(poly, length)
```

Substitute every ``x_k`` in ``poly`` by its ``fgp_operator`` (quantize a polynomial).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.mul_scalar"></a>

#### mul\_scalar

```python
def mul_scalar(elem, other)
```

Multiply on the right by a polynomial/Schubert element, commuting it past the operators via
``schubmult_double_down`` (Leibniz rule for divided differences).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.mul_perm"></a>

#### mul\_perm

```python
def mul_perm(elem, perm)
```

Right-multiply every operator ``partial_k`` by ``partial_perm``, keeping only length-additive products.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Left-multiply by a scalar/polynomial (coefficients sit on the left, so this is plain scaling).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Ring product: scalars scale, nilHecke elements combine via ``mul_scalar`` then ``mul_perm``, else ``mul_scalar``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list (the operator ``partial_w``) or a polynomial (a scalar).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The display symbol ``df(w)`` / ``∂(w)`` / ``\partial^w`` for the operator indexed by ``k``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce ``element`` into the coefficient domain, refusing ring elements and anything containing an ``x`` variable.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.genset"></a>

#### genset

```python
@property
def genset()
```

The ``x`` alphabet.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.from_expr"></a>

#### from\_expr

```python
def from_expr(x)
```

Build the scalar element ``x * partial_id``.

