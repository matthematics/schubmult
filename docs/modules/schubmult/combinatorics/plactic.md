<a id="schubmult.combinatorics.plactic"></a>

# schubmult.combinatorics.plactic

`Plactic`: skew semistandard tableaux with the plactic (Knuth) crystal structure, stored as a
grid with an optional inner (skew) shape of holes.

<a id="schubmult.combinatorics.plactic.Plactic"></a>

## Plactic Objects

```python
class Plactic(GridPrint, CrystalGraph)
```

A (skew) semistandard Young tableau stored as a grid, with the ``gl_n`` plactic crystal
structure (Knuth relations). ``inner_shape`` (a partition of leading holes per row) makes it a
skew tableau; construct directly from a grid or via classmethods like ``Plactic.yamanouchi``.

<a id="schubmult.combinatorics.plactic.Plactic.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(row, col)
```

Perform an upward jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new Plactic tableau.

<a id="schubmult.combinatorics.plactic.Plactic.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(row, col)
```

Perform a jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new Plactic tableau.

<a id="schubmult.combinatorics.plactic.Plactic.shiftup"></a>

#### shiftup

```python
def shiftup(k)
```

Return a new Plactic with all entries increased by k.

<a id="schubmult.combinatorics.plactic.Plactic.all_ss_tableaux"></a>

#### all\_ss\_tableaux

```python
@classmethod
def all_ss_tableaux(cls, shape, max_entry, inner_shape=None)
```

Generate all semistandard tableaux of given shape (or skew shape) with entries <= max_entry.

**Arguments**:

  
- `shape` - Sequence of row lengths (outer shape)
- `max_entry` - Maximum entry value
- `inner_shape` - Optional sequence of left offsets per row (for skew shapes).
  If provided, positions [row][0:inner_shape[row]] are marked as 0.
  

**Returns**:

  Set of Plactic instances representing all valid semistandard tableaux

<a id="schubmult.combinatorics.plactic.Plactic.row_word"></a>

#### row\_word

```python
@property
def row_word()
```

Return the row-reading word as a flat tuple.

<a id="schubmult.combinatorics.plactic.Plactic.transpose"></a>

#### transpose

```python
def transpose()
```

Return the transpose of this Plactic tableau.

<a id="schubmult.combinatorics.plactic.Plactic.invert"></a>

#### invert

```python
def invert()
```

Return a Plactic whose entries are remapped so that standard
(increasing) insertion order applies. If reverse_semistandard is True
we negate entries (so larger original becomes smaller).

<a id="schubmult.combinatorics.plactic.Plactic.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other)
```

Plactic product: insert entries of `other` in row-reading order
(top-to-bottom, left-to-right) into a copy of self.

<a id="schubmult.combinatorics.plactic.Plactic.skew_shape"></a>

#### skew\_shape

```python
@property
def skew_shape()
```

Return the skew shape as a tuple of (row_length, left_offset) pairs.

<a id="schubmult.combinatorics.plactic.Plactic.rs_insert"></a>

#### rs\_insert

```python
def rs_insert(*letters)
```

Insert one or more letters in sequence (row-insertion) and return a new Plactic.

<a id="schubmult.combinatorics.plactic.Plactic.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

Crystal raising operator e_i on the Plactic tableau (delegates to RCGraph).

<a id="schubmult.combinatorics.plactic.Plactic.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i)
```

Crystal lowering operator f_i on the Plactic tableau (delegates to RCGraph).

<a id="schubmult.combinatorics.plactic.Plactic.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Return the crystal weight of this tableau (delegated to RCGraph).

<a id="schubmult.combinatorics.plactic.Plactic.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Return the length/number of rows used for the crystal

<a id="schubmult.combinatorics.plactic.Plactic.yamanouchi"></a>

#### yamanouchi

```python
@classmethod
def yamanouchi(cls, shape)
```

Return the Yamanouchi (highest-weight) tableau of the given shape.

<a id="schubmult.combinatorics.plactic.Plactic.is_increasing"></a>

#### is\_increasing

```python
@property
def is_increasing()
```

Check if the tableau is strictly increasing in rows and columns.

<a id="schubmult.combinatorics.plactic.Plactic.reverse_rsk"></a>

#### reverse\_rsk

```python
def reverse_rsk(recording_tableau)
```

Inverse RSK (row-insertion) for the pair (P,Q) where `self` is P and
`recording_tableau` is the standard recording tableau Q of the same shape.

Returns the original word as a list of integers (in insertion order).

<a id="schubmult.combinatorics.plactic.Plactic.rsk_insert"></a>

#### rsk\_insert

```python
@classmethod
def rsk_insert(cls, *letters)
```

Perform ordinary RSK (row insertion) on the given sequence of letters,
starting from this Plactic as the initial P-tableau. Returns a pair
(P_tableau, Q_tableau) where both are Plactic instances and Q is the
standard recording tableau with entries 1..m (in insertion order).

Usage:
  P, Q = Plactic().rsk_insert(3,1,2,1)
  or
  P, Q = Plactic().rsk_insert([3,1,2,1])

<a id="schubmult.combinatorics.plactic.Plactic.reverse_rectify_to_outer"></a>

#### reverse\_rectify\_to\_outer

```python
def reverse_rectify_to_outer(outer_shape)
```

Deterministic reverse-rectification to a given outer shape `outer_shape`.

Given a (straight) tableau `self` of shape lambda, produce a skew tableau
(represented as a Plactic whose rows may contain 0's for inner cells)
of outer shape `outer_shape` whose rectification is `self`.

outer_shape: iterable of nonnegative ints giving the desired outer row
lengths (mu_0 >= mu_1 >= ...).

Algorithm (deterministic):

- Let mu be the set of cells (r,c) with 0 <= r < len(mu) and 0 <= c < mu[r].
- While the current set of occupied cells (from the working tableau) is
a strict subset of mu:

* choose an outer corner cell (r,c) in mu\current_cells (no cell of mu
to its right or below). Choose the maximal such (r,c) (deterministic).
* create a hole at (r,c) (extend rows/cols as needed, set that cell to 0),
then perform an upward jeu-de-taquin slide from (r,c) using
up_jdt_slide to move the hole inward.
* adopt the resulting tableau and continue.

- Return the resulting Plactic (with zeros marking inner/removed cells).

**Notes**:

  
  - Raises ValueError if outer_shape does not dominate the current shape
  (i.e. mu must contain the current occupied cells).
  - Raises RuntimeError if no suitable outer corner can be found or if an
  up_jdt_slide fails (this indicates the requested outer shape is not
  attainable by reverse-rectification).

