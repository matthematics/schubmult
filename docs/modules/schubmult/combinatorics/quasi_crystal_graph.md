<a id="schubmult.combinatorics.quasi_crystal_graph"></a>

# schubmult.combinatorics.quasi\_crystal\_graph

Quasi-crystal variant of `CrystalGraph`: raising/lowering operators (``quasi_raising_operator``/
``quasi_lowering_operator``) that may be undefined even mid-string, used where the standard
crystal axioms only hold up to the extra ``ep[-2][0] > 0 and ep[-1][1] > 0`` obstruction
checked in `QuasiCrystalGraphTensor`.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph"></a>

## QuasiCrystalGraph Objects

```python
class QuasiCrystalGraph(CrystalGraph)
```

Abstract base for quasi-crystal elements; see the module docstring.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.quasi_raising_operator"></a>

#### quasi\_raising\_operator

```python
def quasi_raising_operator(index)
```

The raising operator for the crystal graph.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.quasi_lowering_operator"></a>

#### quasi\_lowering\_operator

```python
def quasi_lowering_operator(index)
```

The lowering operator for the crystal graph.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.to_quasi_lowest_weight"></a>

#### to\_quasi\_lowest\_weight

```python
def to_quasi_lowest_weight()
```

Return the lowest weight element in the same quasi-crystal component.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.reverse_quasi_lower_seq"></a>

#### reverse\_quasi\_lower\_seq

```python
def reverse_quasi_lower_seq(seq)
```

Apply the reverse of the given lowering sequence.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor"></a>

## QuasiCrystalGraphTensor Objects

```python
class QuasiCrystalGraphTensor(QuasiCrystalGraph)
```

Tensor product of quasi-crystal elements (`factors`); mirrors `CrystalGraphTensor` but the
``quasi_lowering_operator``/``quasi_raising_operator`` calls additionally return ``None`` when
``ep[-2][0] > 0 and ep[-1][1] > 0`` at the last two left-folded ``(epsilon, phi)`` entries.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Sum of the factors' weights (zero-padded to the longest).

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.weight_bump"></a>

#### weight\_bump

```python
def weight_bump()
```

Apply ``weight_bump`` to every factor.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.all_highest_weights"></a>

#### all\_highest\_weights

```python
def all_highest_weights()
```

All highest-weight tensors reachable by taking the highest weight of independently chosen
elements from each factor's full quasi-crystal.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.__init__"></a>

#### \_\_init\_\_

```python
def __init__(*factors)
```

Build the tensor product of the given quasi-crystal elements, left to right.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

The maximum ``crystal_length`` over all factors.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.quasi_lowering_operator"></a>

#### quasi\_lowering\_operator

```python
def quasi_lowering_operator(index)
```

Apply ``quasi_lowering_operator(index)`` to the rightmost eligible factor, or ``None``
if the quasi-crystal obstruction blocks it.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.quasi_raising_operator"></a>

#### quasi\_raising\_operator

```python
def quasi_raising_operator(index)
```

Apply ``quasi_raising_operator(index)`` to the rightmost eligible factor, or ``None``
if the quasi-crystal obstruction blocks it.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.epsilon"></a>

#### epsilon

```python
def epsilon(i)
```

``epsilon_i`` of the tensor, read off the last entry of the left-folded ``(epsilon, phi)`` table.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraphTensor.phi"></a>

#### phi

```python
def phi(i)
```

``phi_i`` of the tensor, read off the last entry of the left-folded ``(epsilon, phi)`` table.

