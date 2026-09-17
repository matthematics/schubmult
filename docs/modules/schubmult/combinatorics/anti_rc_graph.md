<a id="schubmult.combinatorics.anti_rc_graph"></a>

# schubmult.combinatorics.anti\_rc\_graph

`AntiRCGraph`: RC graphs viewed "anti" (rows counted from the bottom, entries at least
their anti row label), used by `RCGraph.left_squash`/`squash_decomp` to peel a Grassmannian
factor off the top of a general RC graph.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph"></a>

## AntiRCGraph Objects

```python
class AntiRCGraph(SchubertMonomialGraph, GridPrint, CrystalGraph)
```

An RC graph in "anti" orientation: row ``i`` (1-indexed from the bottom) holds reflections
``>= i``. ``to_rc_graph``/``from_rc_graph`` convert to/from the ordinary `RCGraph` orientation
(row reversal); most other operations (crystal operators, products, squashing) are defined by
delegating to the `RCGraph` view.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.__init__"></a>

#### \_\_init\_\_

```python
def __init__(rows_or_grid: Iterable[Iterable[int]] | np.ndarray,
             *,
             _is_copy: bool = False) -> None
```

Build from a sequence of rows (each an iterable of reflection labels) or a raw 0/1 grid.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.rows"></a>

#### rows

```python
@property
def rows() -> int
```

Number of rows.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.cols"></a>

#### cols

```python
@property
def cols() -> int
```

Number of columns.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.copy"></a>

#### copy

```python
def copy() -> AntiRCGraph
```

Shallow copy.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.has_element"></a>

#### has\_element

```python
def has_element(i: int, j: int) -> bool
```

Whether the reflection is marked at 1-indexed grid position ``(i, j)``.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.perm"></a>

#### perm

```python
@property
def perm() -> Permutation
```

The permutation induced by this anti RC graph: ``~anti_permutation``.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.reflection_view"></a>

#### reflection\_view

```python
@property
def reflection_view() -> RCGraph
```

The rows reversed into ordinary `RCGraph` orientation.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph() -> RCGraph
```

Alias for ``reflection_view``.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc: RCGraph) -> AntiRCGraph
```

Inverse of ``to_rc_graph``: reverse the rows of an ordinary `RCGraph`.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.from_reduced_anticompatible"></a>

#### from\_reduced\_anticompatible

```python
@classmethod
def from_reduced_anticompatible(cls,
                                word: Sequence[int],
                                seq: Sequence[int],
                                length: int | None = None) -> AntiRCGraph
```

Build from a reduced word and its anti-compatible sequence (dual of `RCGraph.from_reduced_compatible`).

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.normalize"></a>

#### normalize

```python
def normalize() -> AntiRCGraph
```

Drop trailing empty rows (via the `RCGraph` view).

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x, y=None, **_kwargs) -> Expr
```

Monomial (or, with ``y``, double) contribution of this anti RC graph to a Schubert polynomial.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.left_zero_act"></a>

#### left\_zero\_act

```python
def left_zero_act() -> set[AntiRCGraph]
```

Set of anti RC graphs obtained from applying the zero-action to the `RCGraph` view.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.right_zero_act"></a>

#### right\_zero\_act

```python
def right_zero_act() -> set[AntiRCGraph]
```

Alias for ``left_zero_act`` (the anti orientation swaps left/right).

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.antiaut"></a>

#### antiaut

```python
def antiaut() -> AntiRCGraph
```

Reverse the row order (an anti-automorphism of the grid).

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.vertical_cut"></a>

#### vertical\_cut

```python
def vertical_cut(row: int) -> tuple[AntiRCGraph, AntiRCGraph]
```

Split at ``row`` into two anti RC graphs (order swapped relative to `RCGraph.vertical_cut`).

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.product"></a>

#### product

```python
def product(other: SchubertMonomialGraph) -> dict[AntiRCGraph, int]
```

RC graph product of ``other`` (stacked above) and ``self``, converted back to anti orientation.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(row: int) -> AntiRCGraph | None
```

Crystal lowering operator, realized via the raising operator of the `RCGraph` view at the mirrored row.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(row: int) -> AntiRCGraph | None
```

Crystal raising operator, realized via the lowering operator of the `RCGraph` view at the mirrored row.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.max_reflection"></a>

#### max\_reflection

```python
@property
def max_reflection() -> int
```

Largest reflection label appearing in any row.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.disjoint_union"></a>

#### disjoint\_union

```python
def disjoint_union(anti_rc: AntiRCGraph) -> AntiRCGraph
```

Stack ``anti_rc`` above ``self`` (shifted so their reflections don't collide), keeping the
same number of rows.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.squash_product"></a>

#### squash\_product

```python
def squash_product(anti_rc: AntiRCGraph) -> AntiRCGraph
```

Product used by `RCGraph.left_squash`: disjoint-union then cut back down to ``self``'s row count.

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.squash_decomp"></a>

#### squash\_decomp

```python
def squash_decomp()
```

Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass.

