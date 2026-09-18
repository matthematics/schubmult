<a id="schubmult.rings.schubert.quantum_schubert_ring"></a>

# schubmult.rings.schubert.quantum\_schubert\_ring

Quantum (single) Schubert polynomial ring: the ``QSx`` interface.

`QuantumSingleSchubertRing` is a `QuantumDoubleSchubertRing` with a zero
coefficient alphabet. This module also re-exports the quantum double and
parabolic quantum rings (``QDSx``, ``QPSx``, ``QPDSx``) for convenience.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing"></a>

## QuantumSingleSchubertRing Objects

```python
class QuantumSingleSchubertRing(QuantumDoubleSchubertRing)
```

The ring of quantum Schubert polynomials ``S^q_w(x)``; ``QSx`` is the standard instance.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.quantize"></a>

#### quantize

```python
def quantize(poly)
```

Quantize a polynomial: expand in classical Schubert polynomials, reinterpret each ``S_w`` as
the quantum ``S^q_w``, and expand back to a polynomial.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants: ``schubmult_q_fast`` when ``basis2`` is this ring, else ``schubmult_q_double_fast``.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product``.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression: single ``x`` variables via ``mult_poly_q``, ``Add``/``Mul``/``Pow``
recursively, anything else as a coefficient.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, a classical or parabolic element (converted
to the quantum basis), or a polynomial expression.

