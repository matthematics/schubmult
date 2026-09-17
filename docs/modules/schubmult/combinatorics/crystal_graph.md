<a id="schubmult.combinatorics.crystal_graph"></a>

# schubmult.combinatorics.crystal\_graph

Abstract Kashiwara/Demazure crystal graph interface, plus the dual, reversed, and
tensor-product crystal constructions built generically on top of it.

Subclasses (`RCGraph`, `Plactic`, `SetLetter`, ...) implement `raising_operator`,
`lowering_operator`, `crystal_weight`, and `crystal_length`; everything else here
(``epsilon``/``phi``, highest/lowest weight, full crystal enumeration, tensor
products) is generic in terms of those four.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph"></a>

## CrystalGraph Objects

```python
class CrystalGraph(Printable)
```

Abstract base for Kashiwara/Demazure crystal elements.

Concrete subclasses must implement `raising_operator`, `lowering_operator`,
`crystal_weight`, and `crystal_length`; operators return ``None`` when undefined
(not an error). Everything else (``epsilon``/``phi``, highest/lowest weight,
full crystal enumeration) is derived generically from those four.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

The raising operator for the crystal graph.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

The lowering operator for the crystal graph.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

The weight of the crystal graph element.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.phi"></a>

#### phi

```python
def phi(i)
```

``phi_i``: the number of times ``lowering_operator(i)`` can be applied before hitting ``None``.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.epsilon"></a>

#### epsilon

```python
def epsilon(i)
```

``epsilon_i``: the number of times ``raising_operator(i)`` can be applied before hitting ``None``.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight(length=None)
```

Return the lowest weight element in the connected component.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight(length=None, start=1)
```

Return the highest weight element in the connected component.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Return the length of the crystal element.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.square"></a>

#### square

```python
@property
def square()
```

Wrap ``self`` so each raising/lowering operator call applies the underlying operator twice
(the "squared" crystal used to detect index-2 symmetrized structures).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply ``lowering_operator`` along ``raise_seq`` in reverse: undo a recorded raising sequence
(e.g. from ``to_highest_weight``) to recover the original element.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.reverse_lower_seq"></a>

#### reverse\_lower\_seq

```python
def reverse_lower_seq(lower_seq)
```

Apply ``raising_operator`` along ``lower_seq`` in reverse: undo a recorded lowering sequence
(e.g. from ``to_lowest_weight``) to recover the original element.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

The ``sl_2``-string reflection at ``index``: raise or lower ``|epsilon_i - phi_i|`` times to
the opposite end of the ``i``-string through ``self``.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.full_crystal"></a>

#### full\_crystal

```python
@property
def full_crystal()
```

All elements reachable from ``self``'s highest-weight element by lowering operators.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.full_crystal_bothways"></a>

#### full\_crystal\_bothways

```python
def full_crystal_bothways(condition=None)
```

All elements reachable from ``self`` by raising or lowering operators in either direction,
optionally restricted to elements satisfying ``condition``.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.full_squared_crystal"></a>

#### full\_squared\_crystal

```python
@property
def full_squared_crystal()
```

All elements reachable from ``self`` via *pairs* of raising/lowering steps (index-2 substructure).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_beneath"></a>

#### crystal\_beneath

```python
@property
def crystal_beneath()
```

All elements reachable from ``self`` by repeated lowering operators.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_above"></a>

#### crystal\_above

```python
def crystal_above(length=None)
```

All elements reachable from ``self`` by repeated raising operators (indices ``1..length-1``).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.truncated_crystal"></a>

#### truncated\_crystal

```python
def truncated_crystal(length, start=1)
```

All elements reachable from ``self``'s highest-weight element (computed at ``length``) by
lowering operators restricted to indices ``start..length-1``.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.params"></a>

#### params

```python
@property
def params()
```

For each index ``i``, the number of times ``raising_operator(i)`` can be applied
(the ``i``-string parameters of ``self``, a poor man's weight-in-a-box coordinate).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.weight_bump"></a>

#### weight\_bump

```python
def weight_bump()
```

Hook for subclasses: an element with the same crystal structure but perturbed so that
``crystal_reflection`` is guaranteed to succeed (used as a fallback by ``weight_reflection``).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.weight_reflection"></a>

#### weight\_reflection

```python
def weight_reflection(index)
```

Like ``crystal_reflection``, but falls back to ``weight_bump`` first if the direct reflection fails.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.is_highest_weight"></a>

#### is\_highest\_weight

```python
@property
def is_highest_weight()
```

Whether no raising operator is defined on ``self`` (it is highest weight in its component).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.is_lowest_weight"></a>

#### is\_lowest\_weight

```python
@property
def is_lowest_weight()
```

Whether no lowering operator is defined on ``self`` (it is lowest weight in its component).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.dual"></a>

#### dual

```python
@property
def dual()
```

Wrap ``self`` in `CrystalGraphDual` (raising/lowering operators swapped, weight negated).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.reverse"></a>

#### reverse

```python
@property
def reverse()
```

Wrap ``self`` in `CrystalGraphReverse` (indices reversed, ``i <-> length + 1 - i``).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphDual"></a>

## CrystalGraphDual Objects

```python
class CrystalGraphDual(CrystalGraph)
```

Dual crystal: raising and lowering operators are swapped and the weight is negated.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphReverse"></a>

## CrystalGraphReverse Objects

```python
class CrystalGraphReverse(CrystalGraph)
```

Crystal with indices reversed: index ``i`` acts as index ``n + 1 - i`` on the base crystal.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor"></a>

## CrystalGraphTensor Objects

```python
class CrystalGraphTensor(CrystalGraph)
```

Tensor product of crystal elements (`factors`), with the standard (signature-rule) tensor
product crystal structure -- raising/lowering act on the leftmost/rightmost eligible factor
according to the left-folded ``(epsilon, phi)`` values (``_left_folded_ep_phi``).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Sum of the factors' weights (zero-padded to the longest).

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.weight_bump"></a>

#### weight\_bump

```python
def weight_bump()
```

Apply ``weight_bump`` to every factor.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.all_highest_weights"></a>

#### all\_highest\_weights

```python
def all_highest_weights()
```

All highest-weight tensors reachable by taking the highest weight of independently chosen
elements from each factor's full crystal.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.__init__"></a>

#### \_\_init\_\_

```python
def __init__(*factors)
```

Build the tensor product of the given crystal elements, left to right.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

The maximum ``crystal_length`` over all factors.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Apply ``lowering_operator(index)`` to the rightmost factor for which the tensor signature rule allows it.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Apply ``raising_operator(index)`` to the rightmost factor for which the tensor signature rule allows it.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.epsilon"></a>

#### epsilon

```python
def epsilon(i)
```

``epsilon_i`` of the tensor, read off the last entry of the left-folded ``(epsilon, phi)`` table.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraphTensor.phi"></a>

#### phi

```python
def phi(i)
```

``phi_i`` of the tensor, read off the last entry of the left-folded ``(epsilon, phi)`` table.

