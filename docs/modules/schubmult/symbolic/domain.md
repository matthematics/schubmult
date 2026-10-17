<a id="schubmult.symbolic.domain"></a>

# schubmult.symbolic.domain

SymPy-free stand-ins for the parts of SymPy's polys ``Domain`` protocol that `schubmult.rings` uses.

The ring classes only ever relied on ``Domain.__call__`` (construct via ``new``), ``Domain.sum``,
``repr == str``, the ``EXRAW`` coefficient domain's ``zero``/``one``, and ``CoercionFailed``.

<a id="schubmult.symbolic.domain.CoercionFailed"></a>

## CoercionFailed Objects

```python
class CoercionFailed(Exception)
```

Raised when a value cannot be coerced into a ring or its coefficient domain.

<a id="schubmult.symbolic.domain.DomainElement"></a>

## DomainElement Objects

```python
class DomainElement()
```

Marker base class for ring elements.

<a id="schubmult.symbolic.domain.Ring"></a>

## Ring Objects

```python
class Ring()
```

Base class for rings: calling a ring constructs an element via ``new``.

<a id="schubmult.symbolic.domain.CompositeDomain"></a>

## CompositeDomain Objects

```python
class CompositeDomain()
```

Marker base class for rings built over a coefficient domain.

