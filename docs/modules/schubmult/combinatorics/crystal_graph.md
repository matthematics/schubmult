<a id="schubmult.combinatorics.crystal_graph"></a>

# schubmult.combinatorics.crystal\_graph

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph"></a>

## CrystalGraph Objects

```python
class CrystalGraph(Printable)
```

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

