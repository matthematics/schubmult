<a id="schubmult.combinatorics.hpd"></a>

# schubmult.combinatorics.hpd

Bumpless Pipe Dreams (HPD) module

<a id="schubmult.combinatorics.hpd.HPDTile"></a>

## HPDTile Objects

```python
class HPDTile(IntEnum)
```

Enumeration of the possible tile types in a pipe dream.

Each tile represents how two pipes (horizontal and vertical) interact in a square.
Note: Whether a tile is "weighty" depends on the row's _id_vector value, not the tile itself.

<a id="schubmult.combinatorics.hpd.HPDTile.TBD"></a>

#### TBD

Placeholder for uninitialized tile

<a id="schubmult.combinatorics.hpd.HPDTile.BLANK"></a>

#### BLANK

Both pipes go straight (no crossing, no elbow)

<a id="schubmult.combinatorics.hpd.HPDTile.CROSS"></a>

#### CROSS

Pipes cross each other

<a id="schubmult.combinatorics.hpd.HPDTile.HORIZ"></a>

#### HORIZ

Horizontal pipe

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_NW"></a>

#### ELBOW\_NW

Elbow: bottom-right to top-left (╯)

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_SE"></a>

#### ELBOW\_SE

Elbow: top-left to bottom-right (╮)

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_NE"></a>

#### ELBOW\_NE

Elbow: bottom-right to top-left (╯)

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_SW"></a>

#### ELBOW\_SW

Elbow: top-left to bottom-right (╮)

<a id="schubmult.combinatorics.hpd.HPDTile.BUMP"></a>

#### BUMP

Bump/osculating tile (pipes touch at corner)

<a id="schubmult.combinatorics.hpd.HPDTile.__str__"></a>

#### \_\_str\_\_

```python
def __str__() -> str
```

Return base display symbol (use get_display_symbol() for context-aware rendering)

<a id="schubmult.combinatorics.hpd.HPDTile.get_display_symbol"></a>

#### get\_display\_symbol

```python
def get_display_symbol(is_weighty: bool) -> str
```

Get display symbol based on whether the tile is in a weighty row

<a id="schubmult.combinatorics.hpd.HPDTile.from_tiletype"></a>

#### from\_tiletype

```python
@classmethod
def from_tiletype(cls, tile: TileType) -> HPDTile
```

Convert from TileType to HPDTile

<a id="schubmult.combinatorics.hpd.HPDTile.is_crossing"></a>

#### is\_crossing

```python
@cached_property
def is_crossing() -> bool
```

True if this tile is a crossing

<a id="schubmult.combinatorics.hpd.HPDTile.is_elbow"></a>

#### is\_elbow

```python
@cached_property
def is_elbow() -> bool
```

True if this tile is any type of elbow

<a id="schubmult.combinatorics.hpd.HPDTile.is_empty"></a>

#### is\_empty

```python
@cached_property
def is_empty() -> bool
```

True if this tile is empty (pipes go straight)

<a id="schubmult.combinatorics.hpd.HPDTile.feeds_right"></a>

#### feeds\_right

```python
@cached_property
def feeds_right() -> bool
```

True if the horizontal pipe continues to the right

<a id="schubmult.combinatorics.hpd.HPDTile.feeds_up"></a>

#### feeds\_up

```python
@cached_property
def feeds_up() -> bool
```

True if the vertical pipe continues upwards

<a id="schubmult.combinatorics.hpd.HPDTile.entrance_from_bottom"></a>

#### entrance\_from\_bottom

```python
@cached_property
def entrance_from_bottom() -> bool
```

True if a pipe can enter from the bottom

<a id="schubmult.combinatorics.hpd.HPDTile.entrance_from_left"></a>

#### entrance\_from\_left

```python
@cached_property
def entrance_from_left() -> bool
```

True if a pipe can enter from the left

<a id="schubmult.combinatorics.hpd.HPD"></a>

## HPD Objects

```python
class HPD(SchubertMonomialGraph, DefaultPrinting)
```

Bumpless Pipe Dream representation.

A bumpless pipe dream is an n×n grid where:
- HPDTile.CROSS (1) represents a crossing
- HPDTile.BLANK (0) represents an empty box (pipes go straight)
- For general pipe dreams, can use HPDTile.ELBOW_* (2-5) for elbows

Each HPD corresponds to a permutation and has an associated weight.

<a id="schubmult.combinatorics.hpd.HPD.__init__"></a>

#### \_\_init\_\_

```python
def __init__(grid, id_vector, *, _is_copy=False) -> None
```

Initialize a HPD from a grid.

**Arguments**:

- `grid` - n×n array-like of HPDTile values, integers 0-5, or list of lists

<a id="schubmult.combinatorics.hpd.HPD.is_classic_row"></a>

#### is\_classic\_row

```python
def is_classic_row(row: int) -> bool
```

Check if row is a classic row (id_vector[row] == 0)

<a id="schubmult.combinatorics.hpd.HPD.concat"></a>

#### concat

```python
@classmethod
def concat(cls, rc, bpd)
```

Concatenate a BPD and RCGraph into a HPD.

**Arguments**:

- `bpd` - BPD instance
- `rc` - RCGraph instance

<a id="schubmult.combinatorics.hpd.HPD.is_weighty_position"></a>

#### is\_weighty\_position

```python
def is_weighty_position(row: int, col: int) -> bool
```

Check if a position is in; a weighty row (id_vector[row] == 1)

<a id="schubmult.combinatorics.hpd.HPD.row_index_to_label"></a>

#### row\_index\_to\_label

```python
def row_index_to_label(row_index: int) -> int
```

Map physical row index (0-based) to row label.

Row labels are assigned counterclockwise:
- _id_vector == 1 rows: labeled 1, 2, ... on RIGHT side, going UPWARD (bottom to top)
- _id_vector == 0 rows: labeled next, on LEFT side, going DOWNWARD (top to bottom)

**Arguments**:

- `row_index` - Physical row index (0-based, top to bottom)
  

**Returns**:

  Row label (1-based)

<a id="schubmult.combinatorics.hpd.HPD.row_label_to_index"></a>

#### row\_label\_to\_index

```python
def row_label_to_index(label: int) -> int
```

Map row label (1-based) to physical row index (0-based).

Inverse of row_index_to_label.

**Arguments**:

- `label` - Row label (1-based)
  

**Returns**:

  Physical row index (0-based, top to bottom)

<a id="schubmult.combinatorics.hpd.HPD.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc: RCGraph) -> HPD
```

Create a HPD from an RC graph.

**Arguments**:

- `rc` - RCGraph instance

<a id="schubmult.combinatorics.hpd.HPD.swap_rows"></a>

#### swap\_rows

```python
def swap_rows(row: int) -> HPD
```

Swap a classic row with the row below it.

**Arguments**:

- `row` - Index of the classic row to swap (0-based)
  

**Returns**:

  New HPD with the specified rows swapped

<a id="schubmult.combinatorics.hpd.HPD.pipe_source_labels"></a>

#### pipe\_source\_labels

```python
def pipe_source_labels(row: int, col: int) -> dict[str, int | None]
```

Determine which row label(s) the pipe(s) at position (row, col) came from.

Returns a dict with keys 'top', 'bottom', 'left', 'right' indicating which
row label the pipe on each side of the tile belongs to (None if no pipe on that side).

Invariants:
- For non-BUMP/CROSS tiles: exactly 2 non-None sides with equal values
- For BUMP: right == bottom, left == top
- For CROSS: top == bottom, left == right

**Arguments**:

- `row` - Physical row index (0-based)
- `col` - Physical column index (0-based)
  

**Returns**:

  Dict with keys 'top', 'bottom', 'left', 'right' mapping to row labels or None

<a id="schubmult.combinatorics.hpd.HPD.from_bruhat_path"></a>

#### from\_bruhat\_path

```python
@classmethod
def from_bruhat_path(cls, path: Sequence[Permutation]) -> HPD
```

Create a HPD from a Bruhat path.

<a id="schubmult.combinatorics.hpd.HPD.row_from_k_chain"></a>

#### row\_from\_k\_chain

```python
@staticmethod
def row_from_k_chain(u: Permutation, w: Permutation, k: int,
                     n: int) -> np.ndarray
```

Construct a single row of tiles from a k-chain according to Definition 3.15.

Given two permutations u and w where u ≤ w in Bruhat order, finds a maximal
k-chain from u to w and constructs a row of n tiles based on that chain.

**Arguments**:

- `u` - Starting permutation
- `w` - Target permutation (must satisfy u ≤ w in Bruhat order)
- `k` - The chain parameter (k >= 1)
  

**Returns**:

  1D numpy array of HPDTile values representing the row
  
  Definition 3.15 cases (for tile at position (row, c)):
  - If chain swaps c with larger but not smaller: ELBOW_SE (⌜)
  - If chain swaps c with both larger and smaller: CROSS (╋)
  - If chain swaps c with smaller but not larger: ELBOW_NW (⌟)
  - If c not among first k numbers of w: BLANK (□)
  - If chain swaps values a,b with a < c < b: BUMP (╬)
  - Otherwise: CROSS (■)

<a id="schubmult.combinatorics.hpd.HPD.build"></a>

#### build

```python
def build() -> None
```

Build internal structures by resolving TBD tiles using lookup table.

<a id="schubmult.combinatorics.hpd.HPD.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Return the size n of the n×n grid

<a id="schubmult.combinatorics.hpd.HPD.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key) -> HPDTile | np.ndarray
```

Access grid elements, casting to HPDTile

<a id="schubmult.combinatorics.hpd.HPD.shiftup"></a>

#### shiftup

```python
def shiftup(shift: int = 1) -> HPD
```

Shift the HPD up by a given amount.

<a id="schubmult.combinatorics.hpd.HPD.perm"></a>

#### perm

```python
@property
def perm() -> Permutation
```

Compute the permutation associated with this HPD.

The permutation is determined by following each vertical pipe from bottom to top.
Pipes enter from the bottom (vertical) and left (horizontal).

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.hpd.HPD.permutation"></a>

#### permutation

```python
@property
def permutation() -> Permutation
```

Alias for perm property

<a id="schubmult.combinatorics.hpd.HPD.inv"></a>

#### inv

```python
@property
def inv() -> int
```

Return the inversion count of the associated permutation.

This is a convenience property that delegates to perm.inv.

**Returns**:

  Number of inversions in the permutation

<a id="schubmult.combinatorics.hpd.HPD.length_vector"></a>

#### length\_vector

```python
@property
def length_vector() -> tuple[int, ...]
```

Compute the length vector of the permutation represented by this HPD.

The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
of weighty tiles in row i. Which tiles are weighty depends on _id_vector[i]:
- If _id_vector[i] == 1: BLANK tiles are weighty
- If _id_vector[i] == 0: CROSS and HORIZ tiles are weighty

**Returns**:

  Tuple of integers representing the length vector

<a id="schubmult.combinatorics.hpd.HPD.from_asm"></a>

#### from\_asm

```python
@classmethod
def from_asm(cls, asm) -> HPD
```

Create a HPD from an ASM (Alternating Sign Matrix).

**Arguments**:

- `asm` - n×n array-like of integers (-1, 0, 1)

**Returns**:

  HPD object

<a id="schubmult.combinatorics.hpd.HPD.weight"></a>

#### weight

```python
@property
def weight() -> Tuple[int, ...]
```

Compute the weight of this HPD.

The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
of empty squares (0s) in column i.

**Returns**:

  Tuple of integers representing the weight

<a id="schubmult.combinatorics.hpd.HPD.word"></a>

#### word

```python
@property
def word() -> Tuple[int, ...]
```

Compute a reduced word for the permutation represented by this HPD.

For each crossing at position (i,j), the word value is the number of pipes
weakly northeast of the crossing minus 1. Weakly northeast means all positions
(r,c) where r <= i and c >= j.

**Returns**:

  Tuple of integers representing the reduced word (1-indexed positions)

<a id="schubmult.combinatorics.hpd.HPD.set_width"></a>

#### set\_width

```python
def set_width(width)
```

Set the width of the HPD by adding empty columns on the right if needed.

<a id="schubmult.combinatorics.hpd.HPD.is_valid"></a>

#### is\_valid

```python
@property
def is_valid() -> bool
```

Check if this is a valid bpd.

**Returns**:

  True if valid, False otherwise

<a id="schubmult.combinatorics.hpd.HPD.__eq__"></a>

#### \_\_eq\_\_

```python
def __eq__(other: object) -> bool
```

Check equality of two HPDs

<a id="schubmult.combinatorics.hpd.HPD.__hash__"></a>

#### \_\_hash\_\_

```python
def __hash__() -> int
```

Hash for use in sets and dicts

<a id="schubmult.combinatorics.hpd.HPD.copy"></a>

#### copy

```python
def copy() -> HPD
```

Create a copy of this HPD

<a id="schubmult.combinatorics.hpd.HPD.num_crossings"></a>

#### num\_crossings

```python
@property
def num_crossings() -> int
```

Total number of crossings in the HPD

<a id="schubmult.combinatorics.hpd.HPD.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the HPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.hpd.HPD.left_root_at"></a>

#### left\_root\_at

```python
def left_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the HPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.hpd.HPD.monk_insert"></a>

#### monk\_insert

```python
def monk_insert(row)
```

RETURNS NORMALIZED

<a id="schubmult.combinatorics.hpd.HPD.combine"></a>

#### combine

```python
def combine(other, shift=None) -> HPD
```

Shift the HPD up by adding empty rows at the bottom.

<a id="schubmult.combinatorics.hpd.HPD.product"></a>

#### product

```python
def product(other: HPD) -> dict[HPD, int]
```

Compute the product of this HPD with another.

<a id="schubmult.combinatorics.hpd.HPD.prod_with_bpd"></a>

#### prod\_with\_bpd

```python
def prod_with_bpd(other: HPD) -> HPD
```

Deprecated: Use product() instead. Returns the single HPD from product dictionary.

<a id="schubmult.combinatorics.hpd.HPD.rebuild"></a>

#### rebuild

```python
def rebuild() -> None
```

Rebuild the HPD to resolve any TBD tiles

<a id="schubmult.combinatorics.hpd.HPD.snap_width"></a>

#### snap\_width

```python
def snap_width() -> HPD
```

Snap the width of the HPD to the length of its permutation.

<a id="schubmult.combinatorics.hpd.HPD.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x: Sequence[Expr],
              y: Sequence[Expr] | None = None,
              **_kwargs) -> Expr
```

Compute the Schubert polynomial value for this HPD.

**Arguments**:

- `x` - Variable or list of variables for polynomial
- `y` - Optional second set of variables for double Schubert polynomial
- ```**_kwargs``` - Additional keyword arguments for polynomial computation (unused)

<a id="schubmult.combinatorics.hpd.HPD.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph() -> RCGraph
```

Convert this HPD to an RC-graph representation.

**Returns**:

  RCGraph object (if available in the module)

