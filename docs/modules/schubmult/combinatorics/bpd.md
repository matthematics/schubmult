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

<a id="schubmult.combinatorics.bpd.TileType.__str__"></a>

#### \_\_str\_\_

```python
def __str__() -> str
```

Single-character display glyph for the tile.

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

<a id="schubmult.combinatorics.bpd.BPD.as_planar_history"></a>

#### as\_planar\_history

```python
def as_planar_history() -> PlanarHistory
```

Convert to a `PlanarHistory` grid of `Tile` objects (via `TileType.as_tile`).

<a id="schubmult.combinatorics.bpd.BPD.rows"></a>

#### rows

```python
@property
def rows() -> int
```

Number of rows.

<a id="schubmult.combinatorics.bpd.BPD.cols"></a>

#### cols

```python
@property
def cols() -> int
```

Number of columns.

<a id="schubmult.combinatorics.bpd.BPD.clear_unreduced_cache"></a>

#### clear\_unreduced\_cache

```python
@classmethod
def clear_unreduced_cache(cls) -> None
```

Clear the memoization caches used by `all_unreduced_bpds`.

<a id="schubmult.combinatorics.bpd.BPD.all_bpds"></a>

#### all\_bpds

```python
@classmethod
def all_bpds(cls,
             w: Permutation,
             length: int | None = None,
             weight: tuple[int] | None = None) -> set[BPD]
```

All (reduced) BPDs for ``w`` with ``length`` rows, optionally restricted to a given ``weight``.

Built recursively from Bruhat paths (`pull_out_var` chains) via the nested ``bruhat_bpd`` helper;
results are cached by ``(w, length)``/``(w, weight)``.

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

<a id="schubmult.combinatorics.bpd.BPD.delete_top_row"></a>

#### delete\_top\_row

```python
def delete_top_row()
```

Remove the top row, first popping it off with `pop_op` until it's exhausted.

<a id="schubmult.combinatorics.bpd.BPD.append"></a>

#### append

```python
def append(other: BPD) -> BPD
```

Stack ``other`` below ``self`` (analogous to `RCGraph.product` but on BPD grids), tracing
pipes across the boundary to resolve the joining tiles.

<a id="schubmult.combinatorics.bpd.BPD.from_bruhat_path"></a>

#### from\_bruhat\_path

```python
@classmethod
def from_bruhat_path(cls, path: Sequence[Permutation]) -> BPD
```

Create a BPD from a Bruhat path, as per Yu
"Embedding bumpless pipedreams as Bruhat chains" (2024)

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
  
  From Yu "Embedding bumpless pipedreams as Bruhat chains" (2024):
  Definition 3.15 cases (for tile at position (row, c)):
  - If chain swaps c with larger but not smaller: ELBOW_SE (⌜)
  - If chain swaps c with both larger and smaller: CROSS (╋)
  - If chain swaps c with smaller but not larger: ELBOW_NW (⌟)
  - If c not among first k numbers of w: BLANK (□)
  - If chain swaps values a,b with a < c < b: BUMP (╬)
  - Otherwise: CROSS (■)

<a id="schubmult.combinatorics.bpd.BPD.to_bruhat_path"></a>

#### to\_bruhat\_path

```python
def to_bruhat_path()
```

Recover the Bruhat path (one permutation per row cut) that produces this BPD via `from_bruhat_path`.

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

<a id="schubmult.combinatorics.bpd.BPD.co_bpd"></a>

#### co\_bpd

```python
def co_bpd()
```

The complementary BPD: swap HORIZ<->CROSS and VERT<->BLANK, reading rows bottom to top.
See Weigandt, "Changing Bases with Pipe Dream Combinatorics" (2025)

<a id="schubmult.combinatorics.bpd.BPD.groth_to_schub"></a>

#### groth\_to\_schub

```python
@classmethod
@cache
def groth_to_schub(cls, groth_perm: Permutation, beta)
```

Expand the Grothendieck class of ``groth_perm`` in the (beta-deformed) Schubert basis, by
taking the co-BPD of every unreduced BPD and keeping the reduced results.

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
The result is intended to preserve row placement of blanks.

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

Compute the length vector of the BPD.

The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
of blanks (weighty tiles) in row i.

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

<a id="schubmult.combinatorics.bpd.BPD.to_asm"></a>

#### to\_asm

```python
def to_asm()
```

Convert to its alternating sign matrix (``+1``/``-1`` at SE/NW elbows, ``0`` elsewhere).

<a id="schubmult.combinatorics.bpd.BPD.rothe_bpd"></a>

#### rothe\_bpd

```python
@classmethod
@cache
def rothe_bpd(cls, perm: Permutation, num_rows: int | None = None) -> BPD
```

The canonical Rothe BPD for ``perm`` (crossings exactly at the Rothe diagram cells).

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

<a id="schubmult.combinatorics.bpd.BPD.trace_pipe"></a>

#### trace\_pipe

```python
def trace_pipe(i: int, j: int, direction: str | None = None) -> int | None
```

Follow the pipe through cell ``(i, j)`` (entering from ``direction``) out to a grid boundary,
returning the column permutation value it exits at.

<a id="schubmult.combinatorics.bpd.BPD.all_se_elbows"></a>

#### all\_se\_elbows

```python
def all_se_elbows() -> set[tuple[int, int]]
```

All ``(row, col)`` positions holding an SE-elbow tile.

<a id="schubmult.combinatorics.bpd.BPD.all_nw_elbows"></a>

#### all\_nw\_elbows

```python
def all_nw_elbows() -> set[tuple[int, int]]
```

All ``(row, col)`` positions holding an NW-elbow tile.

<a id="schubmult.combinatorics.bpd.BPD.all_blanks"></a>

#### all\_blanks

```python
def all_blanks() -> set[tuple[int, int]]
```

All ``(row, col)`` positions holding a blank tile.

<a id="schubmult.combinatorics.bpd.BPD.all_crossings"></a>

#### all\_crossings

```python
def all_crossings() -> set[tuple[int, int]]
```

All ``(row, col)`` positions holding a crossing tile.

<a id="schubmult.combinatorics.bpd.BPD.all_tiles_of_type"></a>

#### all\_tiles\_of\_type

```python
def all_tiles_of_type(tile_type: TileType) -> set[tuple[int, int]]
```

All ``(row, col)`` positions matching ``tile_type`` (or any type in an iterable of types).

<a id="schubmult.combinatorics.bpd.BPD.droop_moves"></a>

#### droop\_moves

```python
def droop_moves() -> set[tuple[tuple[int, int], tuple[int, int]]]
```

All valid droop moves ``((elbow_pos), (blank_pos))``: legal (SE-elbow, blank) pairs with no
blocking elbow/bump strictly between them.

<a id="schubmult.combinatorics.bpd.BPD.min_droop_moves"></a>

#### min\_droop\_moves

```python
def min_droop_moves() -> set[tuple[tuple[int, int], tuple[int, int]]]
```

The minimal (nearest-blank) droop move available from each SE-elbow/bump, if any.

<a id="schubmult.combinatorics.bpd.BPD.do_min_droop_move"></a>

#### do\_min\_droop\_move

```python
def do_min_droop_move(move: tuple[tuple[int, int], tuple[int, int]]) -> BPD
```

Apply a minimal droop move (from `min_droop_moves`), which may leave a bump tile at the corners.

<a id="schubmult.combinatorics.bpd.BPD.do_droop_move"></a>

#### do\_droop\_move

```python
def do_droop_move(move: tuple[tuple[int, int], tuple[int, int]]) -> BPD
```

Apply a droop move (from `droop_moves`): slide the SE-elbow down-right into the blank corner.

<a id="schubmult.combinatorics.bpd.BPD.huang_bump"></a>

#### huang\_bump

```python
def huang_bump(a, b)
```

Mark the crossing at inversion ``(a, b)`` as a bump and propagate via `_monk_iterate` (BPD
analogue of `RCGraph.huang_bump`).

<a id="schubmult.combinatorics.bpd.BPD.random_bpd"></a>

#### random\_bpd

```python
@classmethod
def random_bpd(cls, perm, num_rows)
```

A uniformly random BPD for ``perm`` with ``num_rows`` rows.

<a id="schubmult.combinatorics.bpd.BPD.monk_insert"></a>

#### monk\_insert

```python
def monk_insert(row)
```

RETURNS NORMALIZED

<a id="schubmult.combinatorics.bpd.BPD.normalize"></a>

#### normalize

```python
def normalize() -> BPD
```

Not yet implemented; intended to trim/pad to the canonical ``len(perm)`` x ``len(perm)`` size.

<a id="schubmult.combinatorics.bpd.BPD.pop_op"></a>

#### pop\_op

```python
def pop_op() -> tuple[BPD, tuple[int, int]]
```

Remove one inversion by the Bergeron-Billey "pop" operation (droop the first blank tile
down-right through the grid); returns ``(new_bpd, (col, row))`` of the popped position.

<a id="schubmult.combinatorics.bpd.BPD.column_perm_at_row"></a>

#### column\_perm\_at\_row

```python
def column_perm_at_row(row: int) -> Permutation
```

Permutation obtained by tracing every pipe entering ``row`` from below out to the bottom boundary.

<a id="schubmult.combinatorics.bpd.BPD.resize"></a>

#### resize

```python
def resize(new_num_rows: int, new_num_cols: int | None = None) -> BPD
```

Grow or shrink to ``new_num_rows`` rows, filling new rows from the Rothe BPD of ``perm``.

<a id="schubmult.combinatorics.bpd.BPD.transpose"></a>

#### transpose

```python
def transpose() -> BPD
```

BPD for ``~perm``: transpose the grid and swap HORIZ/VERT tiles.

<a id="schubmult.combinatorics.bpd.BPD.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph) -> BPD
```

Build the BPD corresponding to an `RCGraph`, via a sequence of `inverse_pop_op` calls
(one per inversion, in reverse reading order).

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

<a id="schubmult.combinatorics.bpd.BPD.inverse_pop_op"></a>

#### inverse\_pop\_op

```python
def inverse_pop_op(*interlaced_rc) -> BPD
```

Inverse of `pop_op`: insert an inversion at ``(col, row)`` pairs, working through each in turn.

<a id="schubmult.combinatorics.bpd.BPD.is_reduced"></a>

#### is\_reduced

```python
@property
def is_reduced()
```

Whether ``self`` is valid and its length vector sums to ``perm.inv`` (no extraneous crossings).

<a id="schubmult.combinatorics.bpd.BPD.as_reduced_compatible"></a>

#### as\_reduced\_compatible

```python
@cache
def as_reduced_compatible()
```

``((col, row), ...)`` pairs recovered by repeatedly applying `pop_op` down to the identity (reversed).

<a id="schubmult.combinatorics.bpd.BPD.rebuild"></a>

#### rebuild

```python
def rebuild() -> None
```

Rebuild the BPD to resolve any TBD tiles

<a id="schubmult.combinatorics.bpd.BPD.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
def zero_out_last_row() -> BPD
```

Drop the last row (`resize` to one fewer row).

<a id="schubmult.combinatorics.bpd.BPD.set_tile"></a>

#### set\_tile

```python
def set_tile(i: int, j: int, tile_type: TileType) -> None
```

Return a copy with the tile at ``(i, j)`` set to ``tile_type``.

<a id="schubmult.combinatorics.bpd.BPD.right_zero_act"></a>

#### right\_zero\_act

```python
def right_zero_act() -> set[BPD]
```

All BPDs one row longer that reduce back to ``self`` under `zero_out_last_row`
(enumerated by trying every subset of new-row crossings and keeping the valid, reduced ones).

<a id="schubmult.combinatorics.bpd.BPD.set_tiles"></a>

#### set\_tiles

```python
def set_tiles(a, b, value: TileType) -> None
```

Return a copy with the tile at ``(a, b)`` set to ``value`` (alias-style variant of ``set_tile``).

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

