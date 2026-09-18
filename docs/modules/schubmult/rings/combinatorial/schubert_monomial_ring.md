<a id="schubmult.rings.combinatorial.schubert_monomial_ring"></a>

# schubmult.rings.combinatorial.schubert\_monomial\_ring

Schubert Monomial Ring module

Provides base classes for rings whose basis elements represent Schubert monomials
(e.g., RC-graphs, BPDs, pipe dreams) with common operations like expansion to
polynomials, divided differences, and crystal operations.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialPrintingTerm"></a>

## SchubertMonomialPrintingTerm Objects

```python
class SchubertMonomialPrintingTerm(TypedPrintingTerm)
```

Printing term for Schubert monomial basis elements.

Delegates printing to the underlying key object (typically an RCGraph, BPD, etc.)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement"></a>

## SchubertMonomialRingElement Objects

```python
class SchubertMonomialRingElement(BaseRingElement)
```

Base class for ring elements whose basis elements are Schubert monomials.

This provides a common interface for objects like:
- RCGraphRingElement (basis elements are RCGraphs)
- BPDRingElement (basis elements are BPDs)

Common operations include:
- Polynomial expansion via polyvalue()
- Divided difference operators
- Crystal structure operations (if the basis elements support them)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x, y=None, *args, **kwargs)
```

Evaluate as a polynomial in variables x (and optionally y).

Linear extension: for each basis element, call its polyvalue() method
and sum the results weighted by coefficients.

**Arguments**:

- `x` - Variable or sequence of variables for polynomial evaluation
- `y` - Optional second set of variables for double Schubert polynomials
- ```**kwargs``` - Additional arguments passed to basis element polyvalue
  

**Returns**:

  Symbolic expression representing the polynomial

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms ``coeff * basis_symbol`` in dict order (sympy printing hook).

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.to_free_algebra_element"></a>

#### to\_free\_algebra\_element

```python
def to_free_algebra_element(basis=None, *, word=False)
```

Convert to FreeAlgebra element in Schubert basis.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing"></a>

## SchubertMonomialRing Objects

```python
class SchubertMonomialRing(BaseRing)
```

Base class for rings whose basis elements are Schubert monomials.

Inherits from BaseRing to provide standard ring operations (add, sub, mul, etc.)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.printing_term"></a>

#### printing\_term

```python
def printing_term(key)
```

Wrap the basis key in a `SchubertMonomialPrintingTerm`.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.from_dict"></a>

#### from\_dict

```python
def from_dict(dct)
```

Build an element from ``{key: coeff}`` without coefficient coercion.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.mul"></a>

#### mul

```python
def mul(a, b)
```

Multiply two elements via each basis key's ``product`` method (which returns ``{key: coeff}``),
or scale by a scalar ``b``.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.rmul"></a>

#### rmul

```python
def rmul(a, b)
```

Scale by the scalar ``b``.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(key)
```

The basis element for ``key`` with coefficient 1.

