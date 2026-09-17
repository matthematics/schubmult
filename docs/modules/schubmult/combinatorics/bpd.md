<a id="schubmult.combinatorics.bpd"></a>

# schubmult.combinatorics.bpd

Bumpless Pipe Dreams (BPD) module

<a id="schubmult.combinatorics.bpd.TileType"></a>

## TileType Objects

```python
class TileType(IntEnum)
```

Enumeration of the 6 possible tile types in a pipe dream.

Each tile represents how two pipes (horizontal and vertical) interact in a square.

<a id="schubmult.combinatorics.bpd.TileType.TBD"></a>

#### TBD

Placeholder for uninitialized tile

<a id="schubmult.combinatorics.bpd.TileType.BLANK"></a>

#### BLANK

Both pipes go straight (no crossing, no elbow)

<a id="schubmult.combinatorics.bpd.TileType.CROSS"></a>

#### CROSS

Pipes cross each other

<a id="schubmult.combinatorics.bpd.TileType.ELBOW_NW"></a>

#### ELBOW\_NW

Elbow: bottom-right to top-left (╯)

<a id="schubmult.combinatorics.bpd.TileType.ELBOW_SE"></a>

#### ELBOW\_SE

Elbow: top-left to bottom-right (╮)

<a id="schubmult.combinatorics.bpd.TileType.BUMP"></a>

#### BUMP

Bump/osculating tile (pipes touch at corner)

<a id="schubmult.combinatorics.bpd.TileType.as_tile"></a>

#### as\_tile

```python
def as_tile() -> Tile
```

Convert TileType to a Tile object with edge information.

<a id="schubmult.combinatorics.bpd.TileType.is_crossing"></a>

#### is\_crossing

```python
@cached_property
def is_crossing() -> bool
```

True if this tile is a crossing

<a id="schubmult.combinatorics.bpd.TileType.is_elbow"></a>

#### is\_elbow

```python
@cached_property
def is_elbow() -> bool
```

True if this tile is any type of elbow

<a id="schubmult.combinatorics.bpd.TileType.is_empty"></a>

#### is\_empty

```python
@cached_property
def is_empty() -> bool
```

True if this tile is empty (pipes go straight)

<a id="schubmult.combinatorics.bpd.TileType.feeds_right"></a>

#### feeds\_right

```python
@cached_property
def feeds_right() -> bool
```

True if the horizontal pipe continues to the right

<a id="schubmult.combinatorics.bpd.TileType.feeds_up"></a>

#### feeds\_up

```python
@cached_property
def feeds_up() -> bool
```

True if the vertical pipe continues upwards

<a id="schubmult.combinatorics.bpd.TileType.entrance_from_bottom"></a>

#### entrance\_from\_bottom

```python
@cached_property
def entrance_from_bottom() -> bool
```

True if a pipe can enter from the bottom

<a id="schubmult.combinatorics.bpd.TileType.entrance_from_left"></a>

#### entrance\_from\_left

```python
@cached_property
def entrance_from_left() -> bool
```

True if a pipe can enter from the left

<a id="schubmult.combinatorics.bpd.BPD"></a>

## BPD Objects

```python
class BPD(SchubertMonomialGraph, DefaultPrinting)
```

Bumpless Pipe Dream representation.

A bumpless pipe dream is an n×n grid where:
- TileType.CROSS (1) represents a crossing
- TileType.BLANK (0) represents an empty box (pipes go straight)
- For general pipe dreams, can use TileType.ELBOW_* (2-5) for elbows

Each BPD corresponds to a permutation and has an associated weight.

<a id="schubmult.combinatorics.bpd.BPD.__init__"></a>

#### \_\_init\_\_

```python
def __init__(grid,
             column_perm: Permutation | None = None,
             *,
             _is_copy=False) -> None
```

Initialize a BPD from a grid.

**Arguments**:

- `grid` - n×n array-like of TileType values, integers 0-5, or list of lists

<a id="schubmult.combinatorics.bpd.BPD.all_unreduced_bpds"></a>

#### all\_unreduced\_bpds

```python
@classmethod
def all_unreduced_bpds(cls,
                       w: Permutation,
                       length: int | None = None,
                       weight: tuple[int] | None = None) -> set[BPD]
```

Enumerate all (possibly unreduced) BPDs for w.

Defers to :meth:`WCGraph.all_wc_graphs`, converts each WCGraph to its
marked bumpless pipe dream, forgets the marks, and reduces to an
ordinary BPD. This is dramatically faster than enumerating ASMs.

<a id="schubmult.combinatorics.bpd.BPD.from_bruhat_path"></a>

#### from\_bruhat\_path

```python
@classmethod
def from_bruhat_path(cls, path: Sequence[Permutation]) -> BPD
```

Create a BPD from a Bruhat path.

<a id="schubmult.combinatorics.bpd.BPD.row_from_k_chain"></a>

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

  1D numpy array of TileType values representing the row
  
  Definition 3.15 cases (for tile at position (row, c)):
  - If chain swaps c with larger but not smaller: ELBOW_SE (⌜)
  - If chain swaps c with both larger and smaller: CROSS (╋)
  - If chain swaps c with smaller but not larger: ELBOW_NW (⌟)
  - If c not among first k numbers of w: BLANK (□)
  - If chain swaps values a,b with a < c < b: BUMP (╬)
  - Otherwise: CROSS (■)

<a id="schubmult.combinatorics.bpd.BPD.build"></a>

#### build

```python
def build() -> None
```

Build internal structures by resolving TBD tiles using lookup table.

<a id="schubmult.combinatorics.bpd.BPD.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Return the size n of the n×n grid

<a id="schubmult.combinatorics.bpd.BPD.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key) -> TileType | np.ndarray
```

Access grid elements, casting to TileType

<a id="schubmult.combinatorics.bpd.BPD.shiftup"></a>

#### shiftup

```python
def shiftup(shift: int = 1) -> BPD
```

Shift the BPD up by a given amount.

<a id="schubmult.combinatorics.bpd.BPD.perm"></a>

#### perm

```python
@property
def perm() -> Permutation
```

Compute the permutation associated with this BPD.

The permutation is determined by following each vertical pipe from bottom to top.
Pipes enter from the bottom (vertical) and left (horizontal).

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.bpd.BPD.permutation"></a>

#### permutation

```python
@property
def permutation() -> Permutation
```

Alias for perm property

<a id="schubmult.combinatorics.bpd.BPD.disjoint_union"></a>

#### disjoint\_union

```python
def disjoint_union(other: BPD) -> BPD
```

Row-preserving disjoint union.

This keeps both summands on the same row indices by placing them side-by-side,
with a single horizontal connector column between them.
The result is intended to preserve row placement of blanks and may be unreduced.

<a id="schubmult.combinatorics.bpd.BPD.disjoint_union_block_diag"></a>

#### disjoint\_union\_block\_diag

```python
def disjoint_union_block_diag(other: BPD) -> BPD
```

Disjoint union via ASM block diagonal (rows of the second summand are shifted).

<a id="schubmult.combinatorics.bpd.BPD.disjoint_union_many"></a>

#### disjoint\_union\_many

```python
@classmethod
def disjoint_union_many(cls, *bpds: BPD) -> BPD
```

Row-preserving disjoint union of any number of BPDs.

<a id="schubmult.combinatorics.bpd.BPD.inv"></a>

#### inv

```python
@property
def inv() -> int
```

Return the inversion count of the associated permutation.

This is a convenience property that delegates to perm.inv.

**Returns**:

  Number of inversions in the permutation

<a id="schubmult.combinatorics.bpd.BPD.length_vector"></a>

#### length\_vector

```python
@property
def length_vector() -> tuple[int, ...]
```

Compute the length vector of the permutation represented by this BPD.

The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
of crossings in row i.

**Returns**:

  Tuple of integers representing the length vector

<a id="schubmult.combinatorics.bpd.BPD.from_asm"></a>

#### from\_asm

```python
@classmethod
def from_asm(cls, asm) -> BPD
```

Create a BPD from an ASM (Alternating Sign Matrix).

**Arguments**:

- `asm` - n×n array-like of integers (-1, 0, 1)

**Returns**:

  BPD object

<a id="schubmult.combinatorics.bpd.BPD.weight"></a>

#### weight

```python
@property
def weight() -> Tuple[int, ...]
```

Compute the weight of this BPD.

The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
of empty squares (0s) in column i.

**Returns**:

  Tuple of integers representing the weight

<a id="schubmult.combinatorics.bpd.BPD.word"></a>

#### word

```python
@property
def word() -> Tuple[int, ...]
```

Compute a reduced word for the permutation represented by this BPD.

For each crossing at position (i,j), the word value is the number of pipes
weakly northeast of the crossing minus 1. Weakly northeast means all positions
(r,c) where r <= i and c >= j.

**Returns**:

  Tuple of integers representing the reduced word (1-indexed positions)

<a id="schubmult.combinatorics.bpd.BPD.set_width"></a>

#### set\_width

```python
def set_width(width)
```

Set the width of the BPD by adding empty columns on the right if needed.

<a id="schubmult.combinatorics.bpd.BPD.is_valid"></a>

#### is\_valid

```python
@property
def is_valid() -> bool
```

Check if this is a valid bpd.

**Returns**:

  True if valid, False otherwise

<a id="schubmult.combinatorics.bpd.BPD.__eq__"></a>

#### \_\_eq\_\_

```python
def __eq__(other: object) -> bool
```

Check equality of two BPDs

<a id="schubmult.combinatorics.bpd.BPD.__hash__"></a>

#### \_\_hash\_\_

```python
def __hash__() -> int
```

Hash for use in sets and dicts

<a id="schubmult.combinatorics.bpd.BPD.copy"></a>

#### copy

```python
def copy() -> BPD
```

Create a copy of this BPD

<a id="schubmult.combinatorics.bpd.BPD.num_crossings"></a>

#### num\_crossings

```python
@property
def num_crossings() -> int
```

Total number of crossings in the BPD

<a id="schubmult.combinatorics.bpd.BPD.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the BPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.bpd.BPD.left_root_at"></a>

#### left\_root\_at

```python
def left_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the BPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.bpd.BPD.monk_insert"></a>

#### monk\_insert

```python
def monk_insert(row)
```

RETURNS NORMALIZED

<a id="schubmult.combinatorics.bpd.BPD.combine"></a>

#### combine

```python
def combine(other, shift=None) -> BPD
```

Shift the BPD up by adding empty rows at the bottom.

<a id="schubmult.combinatorics.bpd.BPD.product"></a>

#### product

```python
def product(other: BPD) -> dict[BPD, int]
```

Compute the product of this BPD with another.

<a id="schubmult.combinatorics.bpd.BPD.prod_with_bpd"></a>

#### prod\_with\_bpd

```python
def prod_with_bpd(other: BPD) -> BPD
```

Deprecated: Use product() instead. Returns the single BPD from product dictionary.

<a id="schubmult.combinatorics.bpd.BPD.rebuild"></a>

#### rebuild

```python
def rebuild() -> None
```

Rebuild the BPD to resolve any TBD tiles

<a id="schubmult.combinatorics.bpd.BPD.snap_width"></a>

#### snap\_width

```python
def snap_width() -> BPD
```

Snap the width of the BPD to the length of its permutation.

<a id="schubmult.combinatorics.bpd.BPD.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x: Sequence[Expr],
              y: Sequence[Expr] | None = None,
              **_kwargs) -> Expr
```

Compute the Schubert polynomial value for this BPD.

**Arguments**:

- `x` - Variable or list of variables for polynomial
- `y` - Optional second set of variables for double Schubert polynomial
- ```**_kwargs``` - Additional keyword arguments for polynomial computation (unused)

<a id="schubmult.combinatorics.bpd.BPD.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph() -> RCGraph
```

Convert this BPD to an RC-graph representation.

**Returns**:

  RCGraph object (if available in the module)

