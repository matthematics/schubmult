<a id="schubmult.rings.combinatorial.plactic_algebra"></a>

# schubmult.rings.combinatorial.plactic\_algebra

`PlacticAlgebra` and `NilPlacticAlgebra`: rings whose basis elements are `Plactic` /
`NilPlactic` tableaux (plactic and nilplactic monoid algebras).

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticPrintingTerm"></a>

## PlacticPrintingTerm Objects

```python
class PlacticPrintingTerm(TypedPrintingTerm)
```

Display symbol for a `PlacticAlgebra` basis tableau.

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticAlgebraElement"></a>

## PlacticAlgebraElement Objects

```python
class PlacticAlgebraElement(BaseRingElement)
```

PlacticAlgebra elements are linear combinations of Plactic basis elements.

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticAlgebra"></a>

## PlacticAlgebra Objects

```python
class PlacticAlgebra(BaseRing)
```

The plactic monoid algebra on `Plactic` tableaux (``op=True`` for the opposite product).

