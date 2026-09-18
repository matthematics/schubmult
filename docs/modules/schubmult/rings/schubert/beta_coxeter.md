<a id="schubmult.rings.schubert.beta_coxeter"></a>

# schubmult.rings.schubert.beta\_coxeter

Scaffold for a beta-deformed (Grothendieck / 0-Hecke) Coxeter operator ring.

Intended to model the beta-isobaric divided differences ``pi_i = partial_i + beta (x_i partial_i - 1)``
and their relations with the simple reflections (see the commented-out relations above
`BetaCoxeterRing`). At present the implementation is an unmodified copy of
`schubmult.rings.schubert.nil_hecke` -- `BetaCoxeterRing`/`BetaCoxeterElement` behave
identically to `NilHeckeRing`/`NilHeckeElement`, and the deformation has not been wired in.
Prefer `nil_hecke` for actual use; this module is kept as a starting point for that work.

<a id="schubmult.rings.schubert.beta_coxeter.BetaCoxeterElement"></a>

## BetaCoxeterElement Objects

```python
class BetaCoxeterElement(DomainElement, DefaultPrinting, dict)
```

An element of a `BetaCoxeterRing`; currently identical in behavior to `NilHeckeElement`.

<a id="schubmult.rings.schubert.beta_coxeter.BetaCoxeterRing"></a>

## BetaCoxeterRing Objects

```python
class BetaCoxeterRing(Ring, CompositeDomain)
```

Scaffold ring; currently identical in behavior to `NilHeckeRing` (see module docstring).

