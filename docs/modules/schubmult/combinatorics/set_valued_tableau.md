<a id="schubmult.combinatorics.set_valued_tableau"></a>

# schubmult.combinatorics.set\_valued\_tableau

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau"></a>

## SetValuedTableau Objects

```python
class SetValuedTableau(GridPrint, CrystalGraph)
```

A semistandard set-valued tableau.

Each box of a Young diagram holds a non-empty *set* of positive integers.
The tableau is semistandard in the set-valued sense: reading each box as its
set of labels, ``max`` of a box is strictly less than ``min`` of the box
below it (columns strictly increase) and weakly less than ``min`` of the box
to its right (rows weakly increase).

Internally the tableau is stored as a dict ``{(row, col): tuple(labels)}``
where each ``labels`` tuple is sorted ascending. Empty boxes are simply
absent from the dict.

The tableau crystal operators are realized via the tensor model on
:class:`~schubmult.combinatorics.set_word.SetWord` with
:class:`~schubmult.combinatorics.set_word.SetLetter` factors, using
row-reading order over boxes (bottom-to-top, left-to-right).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.__init__"></a>

#### \_\_init\_\_

```python
def __init__(cells=None)
```

Create a set-valued tableau.

``cells`` may be:

- a dict ``{(row, col): iterable of labels}``, or
- a nested sequence of rows, where each entry is an iterable of labels
  (or a single integer for a singleton box); ``None`` marks an empty
  leading (skew) box.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.from_cells"></a>

#### from\_cells

```python
@classmethod
def from_cells(cls, cells)
```

Build a :class:`SetValuedTableau` from a dict ``{(row, col): labels}``.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.cells"></a>

#### cells

```python
@property
def cells()
```

Return the underlying ``{(row, col): tuple(labels)}`` dict (a copy).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.weight"></a>

#### weight

```python
@property
def weight()
```

Return the content weight: ``weight[v - 1]`` counts occurrences of ``v``.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.add_label"></a>

#### add\_label

```python
def add_label(row, col, label)
```

Return a new tableau with ``label`` added to the box at ``(row, col)``.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.is_semistandard"></a>

#### is\_semistandard

```python
def is_semistandard()
```

Check the semistandard set-valued conditions on rows and columns.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i)
```

Crystal lowering operator ``f_i`` via :class:`SetWord`.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

Crystal raising operator ``e_i`` via :class:`SetWord`.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Content weight ``(`1`, `2`, ...)`` (alias of :attr:`weight`).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Upper bound on crystal operator indices (matches ``Plactic``).

