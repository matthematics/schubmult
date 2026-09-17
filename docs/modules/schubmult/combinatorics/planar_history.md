<a id="schubmult.combinatorics.planar_history"></a>

# schubmult.combinatorics.planar\_history

Grid of NORTH/EAST/SOUTH/WEST-edged tiles (crossings, bumps, empty cells) tracking pipe
history in a lattice, with conversion to the induced permutation reduced word (``perm_word``).

<a id="schubmult.combinatorics.planar_history.Tile"></a>

## Tile Objects

```python
class Tile()
```

A grid tile, identified by which of its NORTH/EAST/SOUTH/WEST edges are connected (as pairs).

<a id="schubmult.combinatorics.planar_history.PlanarHistory"></a>

## PlanarHistory Objects

```python
class PlanarHistory()
```

A grid of `Tile`s (``CROSS``, ``BUMP``, or ``EMPTY``) recording pipe crossing history;
``perm_word`` reads off the induced permutation word from the crossing positions.

<a id="schubmult.combinatorics.planar_history.PlanarHistory.__init__"></a>

#### \_\_init\_\_

```python
def __init__(grid: np.ndarray)
```

Wrap a 2D array of `Tile` objects.

<a id="schubmult.combinatorics.planar_history.PlanarHistory.rows"></a>

#### rows

```python
@property
def rows()
```

Number of grid rows.

<a id="schubmult.combinatorics.planar_history.PlanarHistory.cols"></a>

#### cols

```python
@property
def cols()
```

Number of grid columns.

<a id="schubmult.combinatorics.planar_history.PlanarHistory.grid"></a>

#### grid

```python
@cached_property
def grid()
```

Copy of the underlying tile array.

<a id="schubmult.combinatorics.planar_history.PlanarHistory.perm_word"></a>

#### perm\_word

```python
@property
def perm_word()
```

The permutation word induced by the ``CROSS`` tiles, read column by column
(bottom-to-top within a column) via the NE pipe-count recurrence.

