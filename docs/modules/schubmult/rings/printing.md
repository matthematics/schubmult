<a id="schubmult.rings.printing"></a>

# schubmult.rings.printing

SymPy atoms used to display ring basis elements.

Every ring's ``printing_term(key)`` returns a `PrintingTerm` subclass instance: an inert SymPy
``Expr`` atom (``args == ()``, so SymPy never traverses into it) that knows how to render itself
for ``str``, pretty printing, and LaTeX. Instances are interned via cached ``__xnew_cached__``
constructors so equal keys give identical objects. The subclasses cover single/double Schubert
(``S``/``DS``), quantum (``QS``/``QDS``, ``QPS``/``QPDS``), Grothendieck (``G``/``DG``), separated
descents (``Xi``), and a `GenericPrintingTerm` ``name(key)`` fallback.

<a id="schubmult.rings.printing.PrintingTerm"></a>

## PrintingTerm Objects

```python
class PrintingTerm(ssymb.Expr)
```

Base display atom carrying a key, generating set, coefficient generating set, and prefix.

<a id="schubmult.rings.printing.GenericPrintingTerm"></a>

## GenericPrintingTerm Objects

```python
class GenericPrintingTerm(PrintingTerm)
```

Displays a key as ``name(key)`` (e.g. ``AGx(perm, n)``, ``N(2, 1)``); the identity key prints as ``1``.

<a id="schubmult.rings.printing.TypedPrintingTerm"></a>

## TypedPrintingTerm Objects

```python
class TypedPrintingTerm(PrintingTerm)
```

Displays a key by delegating to the key's own printer (used for keys that are themselves
printable objects such as RC graphs).

<a id="schubmult.rings.printing.DSchubPoly"></a>

## DSchubPoly Objects

```python
class DSchubPoly(PrintingTerm)
```

Schubert polynomial term: ``S<genset>(perm)`` or ``DS<genset>(perm, <coeff_genset>)``.

<a id="schubmult.rings.printing.SepDescSchubPoly"></a>

## SepDescSchubPoly Objects

```python
class SepDescSchubPoly(PrintingTerm)
```

Separated-descents term for the key ``(perm, numvars)``: ``Xi_{perm}^{numvars}``.

<a id="schubmult.rings.printing.QDSchubPoly"></a>

## QDSchubPoly Objects

```python
class QDSchubPoly(PrintingTerm)
```

Quantum Schubert term: ``QS<genset>(perm)`` or ``QDS<genset>(perm, <coeff_genset>)``.

<a id="schubmult.rings.printing.PQDSchubPoly"></a>

## PQDSchubPoly Objects

```python
class PQDSchubPoly(PrintingTerm)
```

Parabolic quantum Schubert term, tagged with the index composition:
``QPS<genset>(comp)(perm)`` or ``QPDS<genset>(comp)(perm, <coeff_genset>)``.

<a id="schubmult.rings.printing.GrothendieckPoly"></a>

## GrothendieckPoly Objects

```python
class GrothendieckPoly(PrintingTerm)
```

Grothendieck term ``G<genset>(perm)``, or ``G<genset>(perm, numvars)`` when the key carries a variable count.

<a id="schubmult.rings.printing.DoubleGrothendieckPoly"></a>

## DoubleGrothendieckPoly Objects

```python
class DoubleGrothendieckPoly(PrintingTerm)
```

Double Grothendieck term ``DG<genset>(perm, <coeff_genset>)``.

