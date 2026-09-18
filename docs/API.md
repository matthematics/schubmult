<a id="schubmult"></a>

# schubmult

<a id="schubmult.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

Lazy import exported names. This allows `import schubmult` to succeed
even if some optional dependencies are missing.

<a id="schubmult.__dir__"></a>

#### \_\_dir\_\_

```python
def __dir__()
```

Include lazily-exported names in dir()

<a id="schubmult._scripts"></a>

# schubmult.\_scripts

Console and helper scripts shipped with schubmult.

<a id="schubmult._scripts.grothmult_double"></a>

# schubmult.\_scripts.grothmult\_double

<a id="schubmult._scripts.grothmult_double.groth_posify"></a>

#### groth\_posify

```python
def groth_posify(val, var2, var3, msg)
```

Positive FGL form of a coefficient, at ``beta = 1``.

Basis: all products of ``z_i - y_j``, ``(1 + y_j)^{-1}``, ``(1 + z_i)^{-1}``.
Substituting ``u_j = 1 + y_j``, ``v_i = 1 + z_i`` (so ``y = u - 1``,
``z = v - 1``) turns the value into a Laurent polynomial in ``u, v``; clearing
the monomial denominator makes everything polynomial.  Candidates are
difference products times leftover clearing monomials, expanded once, with
``as_coefficients_dict`` terms used as opaque basis vectors for an integer
LP.  ``beta`` is restored as ``beta**(`diffs` - d)`` with the Laurent atoms
``(1 + beta*y)`` degree 0 by construction.

<a id="schubmult._scripts.grothmult_py"></a>

# schubmult.\_scripts.grothmult\_py

<a id="schubmult._scripts.schubmult_double"></a>

# schubmult.\_scripts.schubmult\_double

<a id="schubmult._scripts.schubmult_py"></a>

# schubmult.\_scripts.schubmult\_py

<a id="schubmult._scripts.schubmult_q"></a>

# schubmult.\_scripts.schubmult\_q

<a id="schubmult._scripts.schubmult_q_double"></a>

# schubmult.\_scripts.schubmult\_q\_double

<a id="schubmult.abc"></a>

# schubmult.abc

Ready-made symbols and generating sets, in the spirit of ``sympy.abc``.

Importing this module creates the SymEngine symbols ``x_1..x_99``, ``y_1..y_99``, ``z_1..z_99``,
``q_1..q_99`` and ``β`` and exposes the following names:


- ``x``, ``y``, ``z``, ``q`` -- `GeneratingSet` objects; ``x[i]`` is the symbol ``x_i``. ``x`` is
the default Schubert variable set, ``y``/``z`` the coefficient (double) variable sets, ``q``
the quantum parameters.
- ``beta`` -- the symbol ``β``, the K-theoretic deformation parameter of Grothendieck polynomials.

It also re-exports the unevaluated symmetric polynomial atoms of
`schubmult.symbolic.symmetric_polynomials`:


- ``E`` (alias `FactorialElemSym`) and ``e`` (alias `ElemSym`) -- factorial and ordinary
elementary symmetric polynomials ``E(p, k, xvars, yvars)`` / ``e(p, k, xvars)``.
- ``H`` (alias `FactorialCompleteSym`) and ``h`` (alias `CompleteSym`) -- the complete
homogeneous analogues.
- ``E_q`` (alias `QFactorialElemSym`) -- quantum factorial elementary symmetric polynomials.

```python
from schubmult.abc import x, y, z, q, beta
```
```python
from schubmult.abc import E, e, H, h, E_q
```

**Example**:

  
```python
>>> from schubmult import Sx
>>> from schubmult.abc import x, y, E
>>> Sx([3, 1, 2]).expand()
x_1**2
>>> E(1, 2, [x[1], x[2]], [y[1], y[2]]).expand(func=True)
x_1 + x_2 - y_1 - y_2
```

<a id="schubmult.combinatorics"></a>

# schubmult.combinatorics

Combinatorics package (permutations, RC/BPD/HPD-graphs, tableaux, crystals).

This ``__init__`` currently re-exports nothing; import submodules directly,
e.g. ``from schubmult.combinatorics.permutation import Permutation``. The
commented-out block below is legacy scaffolding from before the module was
split out of ``schub_lib``/``quantum_double``, kept for reference.

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

<a id="schubmult.combinatorics.bpd.BPD.delete_row"></a>

#### delete\_row

```python
def delete_row(row: int) -> BPD
```

Remove ``row`` by tracing its pipe out to the boundary and dropping the corresponding grid row/column.

<a id="schubmult.combinatorics.bpd.BPD.prepend_row"></a>

#### prepend\_row

```python
def prepend_row(value_of_row: int) -> BPD
```

Insert a new top row realizing the given permutation value, growing the permutation by one.

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

<a id="schubmult.combinatorics.chute_move_element"></a>

# schubmult.combinatorics.chute\_move\_element

Chute moves on RC graphs: a marked-row wrapper tracking a before/after pair of RC graphs
related by simultaneous chute moves on a chosen set of non-adjacent rows.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement"></a>

## ChuteMoveElement Objects

```python
class ChuteMoveElement(GridPrint)
```

An RC graph together with the result of applying chute moves at ``rows``.

``rows`` must be pairwise non-adjacent; for each row, moves the element hanging
off the row's end down into the first available gap in the row below,
raising ``ValueError`` if no valid chute move exists. Stores the pair
``(original, moved)`` RC graphs.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.product"></a>

#### product

```python
def product(other)
```

Stack ``self`` above ``other`` (via the underlying RC graph product) and combine their marked rows.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.chute_degree"></a>

#### chute\_degree

```python
@property
def chute_degree()
```

Number of marked rows (simultaneous chute moves applied).

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.chute_move_rows"></a>

#### chute\_move\_rows

```python
@property
def chute_move_rows()
```

The set of marked row indices.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.cols"></a>

#### cols

```python
@property
def cols()
```

Number of columns of the underlying RC graph.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.rows"></a>

#### rows

```python
@property
def rows()
```

Number of rows of the underlying RC graph.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.print_element"></a>

#### print\_element

```python
@property
def print_element()
```

A `GridPrint`-compatible view highlighting the before/after cells at the marked rows.

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

<a id="schubmult.combinatorics.double_forest"></a>

# schubmult.combinatorics.double\_forest

Double forest polynomial construction via the vine subword model.

This module provides a non-script home for the core computational routine
`double_forest_polynomial`, suitable for use by library code.

<a id="schubmult.combinatorics.double_forest.forest_from_code"></a>

#### forest\_from\_code

```python
def forest_from_code(code)
```

Build indexed forest F with c(F)=code via Thompson monoid factorization.

<a id="schubmult.combinatorics.double_forest.forest_qdes"></a>

#### forest\_qdes

```python
def forest_qdes(code)
```

Return the left terminal set (qdes) of the forest with the given code.

qdes(F) = {i : c_i > 0 and c_{i+1} = 0}, using 1-based indexing.

These are the descents of the forest, analogous to descent set of a
permutation (Nadeau-Spink-Tewari, arXiv:2406.01510, §3.2).

<a id="schubmult.combinatorics.double_forest.forest_code_from_trimming_sequence"></a>

#### forest\_code\_from\_trimming\_sequence

```python
def forest_code_from_trimming_sequence(trimming_seq)
```

Return the forest composition sfc(F) for the unique indexed forest F
whose set of trimming sequences contains ``trimming_seq``.

A trimming sequence (i_1, ..., i_k) for F ∈ IndexedForests is defined
recursively: i_k ∈ qdes(F) and (i_1,...,i_{k-1}) ∈ Trim(F/i_k).
(Nadeau-Spink-Tewari, arXiv:2406.01510, Definition 3.8.)

Recovery uses the inverse operation (blossoming):
    sfc(F · i) = (c_1,...,c_{i-1}, c_i+1, 0, c_{i+1}, c_{i+2},...)
i.e. increment position i and insert a 0 immediately after.
Starting from the empty forest and blossoming left-to-right through
(i_1,...,i_k) reconstructs F.

Parameters
----------
trimming_seq : iterable of int
    Sequence of positive integers (1-based leaf positions).

Returns
-------
tuple of int
    The composition sfc(F); trailing zeros are kept to preserve the
    ambient length implied by the sequence.

Examples
--------
>>> forest_code_from_trimming_sequence([1, 1, 2, 4, 7])
(2, 1, 0, 1, 0, 0, 1, 0)
>>> forest_code_from_trimming_sequence([3])
(0, 0, 1, 0)
>>> forest_code_from_trimming_sequence([])
()

<a id="schubmult.combinatorics.double_forest.sylvester_word"></a>

#### sylvester\_word

```python
def sylvester_word(forest)
```

Return one Sylvester word of `forest` via pre-order traversal.

<a id="schubmult.combinatorics.double_forest.sylvester_forest"></a>

#### sylvester\_forest

```python
def sylvester_forest(code, genset, t)
```

Sum of ``polyvalue(genset, t)`` over all RC graphs of ``uncode(code)`` with the given forest weight;
the single (non-double) specialization of `double_sylvester_forest`.

<a id="schubmult.combinatorics.double_forest.double_sylvester_forest"></a>

#### double\_sylvester\_forest

```python
def double_sylvester_forest(code, genset, t)
```

Double (equivariant) forest polynomial, computed by pairing RC graphs of ``u`` and ``v`` from
the double Schubert expansion ``Sx([]) * DSx(perm, "t")`` whose merged vine diagram matches the
principal RC graph's omega-invariant target.

<a id="schubmult.combinatorics.double_forest.canonical_forest_from_word"></a>

#### canonical\_forest\_from\_word

```python
def canonical_forest_from_word(word)
```

Canonical forest form used to test Sylvester-equivalence of words.

<a id="schubmult.combinatorics.double_forest.long_word"></a>

#### long\_word

```python
def long_word(n)
```

Build omega_tilde_[n] as a list of letter records.

<a id="schubmult.combinatorics.double_forest.letter_weight"></a>

#### letter\_weight

```python
def letter_weight(letter, x_gen, t_gen)
```

Weight map in the vine model. x_gen, t_gen are 1-indexed subscriptables.

<a id="schubmult.combinatorics.double_forest.double_forest_polynomial"></a>

#### double\_forest\_polynomial

```python
def double_forest_polynomial(code, x_gen, t_gen, n=None)
```

Compute P_F(x;t) for indexed forest F with code c(F)=code.

<a id="schubmult.combinatorics.double_forest.reflection_subword_polynomial"></a>

#### reflection\_subword\_polynomial

```python
def reflection_subword_polynomial(code, x_gen, t_gen, n=None, mode="forest")
```

Vine subword model from arXiv:2504.15234 (Bergeron-Gagnon-Nadeau-Spink-Tewari).

Iterates subwords pi of the long word `long_word(n)` of size
    sz = sum(code) = ell(uncode(code))
and sums wt(pi) under one of two filters on the *value sequence*
(treating barred and unbarred letters by their value only):

  mode='forest'   :  values must be Sylvester-equivalent to sylvester_word(F),
                     i.e. lie in Syl(F).  Equals P_F(x;t)  (Theorem 5.1).
  mode='schubert' :  values must be a reduced word for perm=uncode(code),
                     i.e. lie in Red(perm).  Equals S_perm(x;t) (Theorem 6.1).

Weights (Sylvester column convention, paper p.19):
    wt(j^(k))      = x_k - t_j
    wt(barred j^(k)) = t_j - t_k
realised by `letter_weight` with bracket-indexed `x_gen`, `t_gen`.

<a id="schubmult.combinatorics.double_forest.reflection_forest_polynomial"></a>

#### reflection\_forest\_polynomial

```python
def reflection_forest_polynomial(code, x_gen, t_gen, n=None, weight_rule=None)
```

Reflection-style forest model from the FULL long word (barred + unbarred).

Iterates subwords of `long_word(n)` (same enumeration used by
`double_forest_polynomial`). A subword is kept iff:
  (a) its value sequence is a reduced expression for the same
      permutation as `sylvester_word(F)`, AND
  (b) its omega-insertion P-symbol (Nadeau-Tewari arXiv:2306.10939 §5.1)
      equals that of `sylvester_word(F)` (single Omega-class).

Default weight is `letter_weight` (paper convention). Override with
`weight_rule(letter, position_in_subword, picked_letters, x_gen, t_gen,
             sylvester_letters)` -> sympy expression, where
    letter            = picked long_word letter dict
    position_in_subword = 0-indexed position in the picked subword
    picked_letters    = tuple of all picked long_word letter dicts
    sylvester_letters = tuple of long_word letters of the canonical
                        Sylvester subword that gave `target`

<a id="schubmult.combinatorics.double_forest.debug_compare_models"></a>

#### debug\_compare\_models

```python
def debug_compare_models(code, x_gen, t_gen, n=None)
```

Print, for `code`, the subwords that survive in:

(A) the alphabet vine model (double_forest_polynomial), and
(B) the reflection-alphabet model with the omega-insertion filter
    using the principal-RC reversed perm-word as target,

so we can stare at the symmetric difference.

<a id="schubmult.combinatorics.hecke_plactic"></a>

# schubmult.combinatorics.hecke\_plactic

`HeckePlactic`: size-preserving Hecke column insertion, a `Plactic` variant.

<a id="schubmult.combinatorics.hecke_plactic.HeckePlactic"></a>

## HeckePlactic Objects

```python
class HeckePlactic(Plactic)
```

Insertion tableau for *size-preserving Hecke column insertion*.

This is a "made up" variant of Hecke (K-theoretic) column insertion that
behaves like Edelman--Greene insertion in that every inserted letter adds
exactly **one** box (the shape/size is preserved: ```boxes` == word length``),
while still preserving the Hecke product / K-Knuth equivalence class of the
reading word.

Ordinary Hecke insertion is *not* size preserving: when a letter cannot
extend the shape it is *absorbed* (no box is added). Here we never absorb --
instead the letter is carried forward to the next column until it can be
placed as a new box. The price is that the resulting tableau ``P`` is only
**column-strict semistandard** (entries strictly increase down columns and
weakly increase along rows) rather than a strict increasing tableau: a value
may repeat within a row.

The recording tableau ``Q`` is an ordinary (single-valued) semistandard
:class:`~schubmult.combinatorics.plactic.Plactic` tableau.

<a id="schubmult.combinatorics.hecke_plactic.HeckePlactic.perm"></a>

#### perm

```python
@property
def perm()
```

Hecke product of the row-reading word (the preserved K-Knuth class).

<a id="schubmult.combinatorics.hecke_plactic.HeckePlactic.hecke_insert"></a>

#### hecke\_insert

```python
def hecke_insert(*letters)
```

Size-preserving Hecke column-insert ``letters`` into a copy of self.

Returns the resulting :class:`HeckePlactic` (the insertion tableau ``P``).
Use :meth:`hecke_insert_rsk` when the recording tableau ``Q`` is required.

<a id="schubmult.combinatorics.hecke_plactic.HeckePlactic.hecke_insert_rsk"></a>

#### hecke\_insert\_rsk

```python
@classmethod
def hecke_insert_rsk(cls, recording, insertion)
```

Size-preserving Hecke column insertion of a two-line array.

``recording`` holds the top-row labels ``k`` (weakly increasing) and
``insertion`` holds the bottom-row letters ``a`` that are column-inserted.

Returns ``(P, Q)`` where ``P`` is a :class:`HeckePlactic` insertion
tableau and ``Q`` is a single-valued semistandard
:class:`~schubmult.combinatorics.plactic.Plactic` recording tableau.
Because insertion is size preserving, ``P`` and ``Q`` have the same shape
and ``Q`` records the label ``k`` in each newly created box.

<a id="schubmult.combinatorics.hecke_plactic.HeckePlactic.reverse_hecke_insert"></a>

#### reverse\_hecke\_insert

```python
def reverse_hecke_insert(corner)
```

Reverse-insert the box at ``corner`` ``(row, col)``.

Returns ``(Y, x)`` where ``Y`` is the resulting :class:`HeckePlactic`
and ``x`` is the reconstructed letter.

<a id="schubmult.combinatorics.hecke_plactic.HeckePlactic.hecke_uninsert_rsk"></a>

#### hecke\_uninsert\_rsk

```python
@classmethod
def hecke_uninsert_rsk(cls, P, Q)
```

Invert :meth:`hecke_insert_rsk`.

Given an insertion tableau ``P`` (:class:`HeckePlactic`) and a recording
tableau ``Q`` (:class:`~schubmult.combinatorics.plactic.Plactic`) of the
same shape, reconstruct the two-line array as ``(recording, insertion)``.

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

<a id="schubmult.combinatorics.increasing_tableau"></a>

# schubmult.combinatorics.increasing\_tableau

`IncreasingTableau`: K-theoretic increasing tableaux (a `Plactic` variant where insertion may
bump without adding a box), plus grid-shape helper utilities.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau"></a>

## IncreasingTableau Objects

```python
class IncreasingTableau(Plactic)
```

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_insert"></a>

#### hecke\_insert

```python
def hecke_insert(*letters)
```

Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.ed_insert_rsk"></a>

#### ed\_insert\_rsk

```python
@classmethod
def ed_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.ed_column_insert_rsk"></a>

#### ed\_column\_insert\_rsk

```python
@classmethod
def ed_column_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(*corners)
```

K-theoretic (Buch–Samuel) jeu de taquin slide toward the top-left.

Starting from one or more empty inner cells, slide entries up/left into
the holes. This is the genuine K-theoretic slide, so a single entry may
migrate into several holes at once; passing multiple corners performs
the simultaneous slide of all of them.

Corners may be given either as ``up_jdt_slide(row, col)`` (a single
corner) or as ``up_jdt_slide((r1, c1), (r2, c2), ...)``. Returns a new
:class:`IncreasingTableau`.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(*corners)
```

K-theoretic (Buch–Samuel) jeu de taquin slide toward the bottom-right.

Starting from one or more empty outer cells, slide entries down/right
into the holes. As with :meth:`up_jdt_slide`, this is the genuine
K-theoretic slide (an entry may migrate into several holes) and passing
multiple corners performs the simultaneous slide of all of them.

Corners may be given either as ``down_jdt_slide(row, col)`` (a single
corner) or as ``down_jdt_slide((r1, c1), (r2, c2), ...)``. Returns a new
:class:`IncreasingTableau`.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.down_jdt_slide_all_inner_corners"></a>

#### down\_jdt\_slide\_all\_inner\_corners

```python
def down_jdt_slide_all_inner_corners()
```

Simultaneously down-slide every valid inner corner.

Collects all inner corners of the skew shape (see
:attr:`iter_inner_corners`) and performs a single K-theoretic down slide
that moves all of them at once. Returns a new
:class:`IncreasingTableau`. If there are no inner corners the tableau is
returned unchanged.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_column_insert"></a>

#### hecke\_column\_insert

```python
def hecke_column_insert(*letters)
```

Column Hecke-insert ``letters`` into a copy of this tableau.

Returns the resulting :class:`IncreasingTableau` (the insertion tableau
``P``). Use :meth:`hecke_column_insert_rsk` when the set-valued
recording tableau ``Q`` is also required.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_column_insert_rsk"></a>

#### hecke\_column\_insert\_rsk

```python
@classmethod
def hecke_column_insert_rsk(cls, recording, insertion)
```

Column Hecke (K-theoretic) insertion of a two-line array.

The two-line array is given by two equal-length sequences: ``recording``
holds the top-row labels ``k`` (weakly increasing, with the letters in
each equal-label block strictly increasing) and ``insertion`` holds the
bottom-row letters ``a`` that are column-inserted.

Returns ``(P, Q)`` where ``P`` is an :class:`IncreasingTableau` (the
insertion tableau) and ``Q`` is a
:class:`~schubmult.combinatorics.set_valued_tableau.SetValuedTableau`,
the semistandard *set-valued* recording tableau.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.reverse_hecke_column_insert"></a>

#### reverse\_hecke\_column\_insert

```python
def reverse_hecke_column_insert(corner, alpha=1)
```

Reverse Hecke-column-insert the box at ``corner`` ``(row, col)``.

``alpha`` is ``1`` when the corner box was created by the forward
insertion (a ``"grow"`` step) and ``0`` when the forward step merely
absorbed a letter at that corner. Returns ``(Y, x)`` where ``Y`` is the
resulting :class:`IncreasingTableau` and ``x`` the reconstructed letter.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_column_uninsert_rsk"></a>

#### hecke\_column\_uninsert\_rsk

```python
@classmethod
def hecke_column_uninsert_rsk(cls, P, Q)
```

Invert :meth:`hecke_column_insert_rsk`.

Given an insertion tableau ``P`` (:class:`IncreasingTableau`) and a
set-valued recording tableau ``Q`` of the same shape, reconstruct the
two-line array as ``(recording, insertion)``.

The recording labels are undone from largest to smallest. Within a box,
the maximum label was the last placed there: a singleton box came from a
``"grow"`` step (reverse with ``alpha = 1``), while a box holding several
labels had its largest label appended by an ``"absorb"`` step (reverse
with ``alpha = 0``, keeping the box). Because insertion is by *columns*
with the canonical convention that equal-label letters are inserted in
strictly *decreasing* order, the corners sharing the maximal label are
undone right-most column first (that being the most recently created
box for the batch).

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other)
```

Plactic product: insert entries of `other` in row-reading order
(top-to-bottom, left-to-right) into a copy of self.

<a id="schubmult.combinatorics.indexed_forests"></a>

# schubmult.combinatorics.indexed\_forests

Indexed forests: binary-search-tree forests indexed by a composition (via the Thompson
monoid factorization), used for the forest basis / forest Schubert-polynomial model
(the ``forest_*`` research scripts in ``_lscripts`` and `schubmult.rings.free_algebra.forest_basis`).

Includes `Node`/`IndexedForest` (the forests themselves), `ParallelInjLetter`/`letterpair`
(parallel-injection alphabet used by omega-insertion), and `LabeledForest`/`DecLabeling`/`LBS`
(labelings of a forest's nodes).

<a id="schubmult.combinatorics.indexed_forests.Node"></a>

## Node Objects

```python
class Node()
```

A node of an indexed forest's binary search tree: an ``index`` (BST key), optional
``label``, and ``left``/``right`` children.

<a id="schubmult.combinatorics.indexed_forests.Node.rho"></a>

#### rho

```python
@property
def rho()
```

Index of the leftmost node in this subtree (its minimum).

<a id="schubmult.combinatorics.indexed_forests.Node.inorder_traversal"></a>

#### inorder\_traversal

```python
@property
def inorder_traversal()
```

Nodes of the subtree in increasing index order.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest"></a>

## IndexedForest Objects

```python
class IndexedForest()
```

A forest of `Node` binary search trees, indexed by a composition (``code``) via the
Thompson monoid factorization (see `forest_from_code`/`double_forest.forest_from_code`).

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.roots"></a>

#### roots

```python
@property
def roots()
```

The root nodes, sorted.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.inorder_traversal"></a>

#### inorder\_traversal

```python
@property
def inorder_traversal()
```

All nodes of the forest in increasing index order.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.code"></a>

#### code

```python
@property
def code()
```

The weak composition indexing this forest (cached).

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.node"></a>

#### node

```python
def node(index)
```

The node with the given index, or ``None``.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.support"></a>

#### support

```python
@property
def support()
```

Sorted tuple of all node indices.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.intervals"></a>

#### intervals

```python
@property
def intervals()
```

Per-root intervals of consecutive indices; see `_forest_intervals`.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.terminal_nodes"></a>

#### terminal\_nodes

```python
@property
def terminal_nodes()
```

Terminal nodes for trimming, indexed by qdes positions.

In Nadeau--Spink--Tewari notation, qdes(F) is defined from the forest
code c(F) by
    qdes(F) = { i : c_i > 0 and c_{i+1} = 0 }.
We realize these as nodes of index i when present in the support.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.trim_descent_nodes"></a>

#### trim\_descent\_nodes

```python
@property
def trim_descent_nodes()
```

Nodes corresponding to trim descents (qdes positions).

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.trim_descents"></a>

#### trim\_descents

```python
@property
def trim_descents()
```

Trim descents (left terminal set qdes) of the indexed forest.

This matches the code-level criterion used in Nadeau--Spink--Tewari:
    qdes(F) = { i : c_i > 0 and c_{i+1} = 0 },
with 1-based indexing and c_{k}=0 beyond the code length.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.trim_descent"></a>

#### trim\_descent

```python
def trim_descent(index)
```

Trim the forest directly at descent ``index``.

The index must lie in ``self.trim_descents``. If ``c = self.code`` and
``index = i`` (1-based), this performs the inverse of blossoming at i:
decrement ``c_i`` by 1 and delete position ``i+1`` (which is 0 for a
valid trim descent), then rebuild the indexed forest from the resulting
code.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.is_left_child"></a>

#### is\_left\_child

```python
def is_left_child(descent)
```

Return whether the leaf labeled ``descent`` is a left child.

The label is interpreted as a support index in the indexed forest.
If the labeled node is not a leaf, this returns ``False``.

<a id="schubmult.combinatorics.indexed_forests.draw_indexed_forest"></a>

#### draw\_indexed\_forest

```python
def draw_indexed_forest(forest,
                        save_path=None,
                        show=True,
                        ax=None,
                        dpi=200,
                        node_size=900,
                        support_color="#1f77b4",
                        edge_color="#333333")
```

Draw an indexed forest with support labels.

Parameters
----------
forest : IndexedForest | iterable[Node]
    Forest to draw. This can be the output of `weak_composition_to_indfor`.
save_path : str | None
    If provided, save the figure (e.g. to a `.png` path).
show : bool
    If True and no external axis is provided, display the figure window.
ax : matplotlib.axes.Axes | None
    Optional existing axis to draw into.
dpi : int
    DPI used when creating/saving a new figure.
node_size : int
    Marker size for node circles.
support_color : str
    Color used for support labels.
edge_color : str
    Color used for edges.

Returns
-------
tuple
    `(fig, ax)` for the rendered drawing.

<a id="schubmult.combinatorics.indexed_forests.ParallelInjLetter"></a>

## ParallelInjLetter Objects

```python
@dataclass(frozen=True, order=True)
class ParallelInjLetter()
```

A letter ``primary[secondary]`` in the parallel-injection alphabet used by omega-insertion.

<a id="schubmult.combinatorics.indexed_forests.make_parallel_injective_word"></a>

#### make\_parallel\_injective\_word

```python
def make_parallel_injective_word(primary_word, secondary_word)
```

Zip two equal-length words into a word of `ParallelInjLetter` biletters.

<a id="schubmult.combinatorics.indexed_forests.weak_composition_to_indfor"></a>

#### weak\_composition\_to\_indfor

```python
def weak_composition_to_indfor(c)
```

Converts a weak composition c+ into an indexed forest.
Returns a list of root nodes for the trees in the forest.

<a id="schubmult.combinatorics.indexed_forests.indexed_forest_from_trimming_word"></a>

#### indexed\_forest\_from\_trimming\_word

```python
def indexed_forest_from_trimming_word(trimming_word)
```

Build an indexed forest directly from a trimming word.

If the trimming word is ``(i_1, ..., i_k)``, we reconstruct the forest code
by iterating the inverse blossoming operation

    c <- (c_1, ..., c_{i-1}, c_i + 1, 0, c_{i+1}, ...)

at each index ``i`` from left to right, then convert the resulting weak
composition code into an ``IndexedForest``.

<a id="schubmult.combinatorics.indexed_forests.minimum_n_for_dual_forest"></a>

#### minimum\_n\_for\_dual\_forest

```python
def minimum_n_for_dual_forest(forest_or_code)
```

Smallest n for which the unsigned ~P_F (eq. 4.1) is nonzero.

With L = internal(F): kappa(v) > rho(v), left child < parent (strict),
right child <= parent (weak); we need labels to fit inside [1..n].

<a id="schubmult.combinatorics.indexed_forests.tilde_forest_polynomial"></a>

#### tilde\_forest\_polynomial

```python
def tilde_forest_polynomial(forest_or_code, n=None)
```

Equation (4.1) of arXiv:2306.10939: the unsigned dual forest polynomial.

    ~P_F = sum_{internal(F)-compatible kappa} x^kappa

Returned as an element of FreeAlgebra(WordBasis), keyed by the
composition (exponent tuple) of length n.

<a id="schubmult.combinatorics.indexed_forests.dual_forest_polynomial"></a>

#### dual\_forest\_polynomial

```python
def dual_forest_polynomial(forest_or_code, n)
```

Signed lower-ideal expansion equal to P_F (Thm. 4.1, arXiv:2306.10939).

    P_F = sum_{lower ideals L of F} (-1)^|L| * sum_{L-compatible kappa} x^kappa

Returned as an element of FreeAlgebra(WordBasis). For comparison against
a directly-computed P_F, use the same n and check equality of the
resulting WordBasis dicts.

<a id="schubmult.combinatorics.indexed_forests.build_balanced_tree"></a>

#### build\_balanced\_tree

```python
def build_balanced_tree(labels)
```

Helper to build a tree where the in-order traversal matches the labels.
This creates the 'canonical labeling' referenced in the paper.

<a id="schubmult.combinatorics.indexed_forests.decreasing_labelings"></a>

#### decreasing\_labelings

```python
def decreasing_labelings(root, max_val, used_vals=None)
```

All labelings of the tree at ``root`` by distinct values ``<= max_val`` that strictly decrease
from each node to its children.

<a id="schubmult.combinatorics.indexed_forests.word_from_labeling"></a>

#### word\_from\_labeling

```python
def word_from_labeling(root, labeling)
```

Word of ``rho`` values of the tree's nodes, ordered by the inverse of the labeling permutation.

<a id="schubmult.combinatorics.indexed_forests.letterpair"></a>

## letterpair Objects

```python
class letterpair()
```

A letter ``primary[secondary]`` (plain-object counterpart of `ParallelInjLetter`, used
where a frozen dataclass isn't convenient).

<a id="schubmult.combinatorics.indexed_forests.LabeledForest"></a>

## LabeledForest Objects

```python
class LabeledForest()
```

An `IndexedForest` together with a labeling of its nodes (``self(index)`` looks up a label).

<a id="schubmult.combinatorics.indexed_forests.DecLabeling"></a>

## DecLabeling Objects

```python
class DecLabeling(LabeledForest)
```

A `LabeledForest` whose labels strictly decrease from each node to its children (checked by ``is_valid``).

<a id="schubmult.combinatorics.indexed_forests.LBS"></a>

## LBS Objects

```python
class LBS(LabeledForest)
```

A `LabeledForest` variant used for the LBS (labeled binary search tree) construction.

<a id="schubmult.combinatorics.indexed_forests.LBS.rootlist"></a>

#### rootlist

```python
@property
def rootlist()
```

The finite part of the rootlist: the root labels together with the bare
letters in the gaps of the support.

The genuine rootlist of Nadeau--Tewari Definition 5.6 also contains every
bare letter below and above the support, so it is infinite in both
directions; use :meth:`rootlist_window` to see those tails.

<a id="schubmult.combinatorics.indexed_forests.LBS.rootlist_window"></a>

#### rootlist\_window

```python
def rootlist_window(lo=None, hi=None)
```

Rootlist of Nadeau--Tewari Definition 5.6, truncated to a window.

Consists of the labels of the roots of the trees -- one for each maximal
interval of the support, not one for each node -- together with the bare
letters ``i[0]`` for which neither ``i`` nor ``i-1`` lies in the support.
The bare part is cofinite, so the genuine rootlist is infinite in both
directions; ``lo`` and ``hi`` bound the range of bare letters reported,
defaulting to one step beyond the support on either side.

<a id="schubmult.combinatorics.indexed_forests.LBS.separators"></a>

#### separators

```python
def separators(a, b, lo=None, hi=None)
```

Rootlist elements lying strictly between the letters ``a`` and ``b``.

<a id="schubmult.combinatorics.indexed_forests.LBS.is_separated"></a>

#### is\_separated

```python
def is_separated(a, b)
```

The criterion of Nadeau--Tewari Proposition 5.8 for ``ab <-> ba``.

<a id="schubmult.combinatorics.indexed_forests.word_to_pair_labeled"></a>

#### word\_to\_pair\_labeled

```python
def word_to_pair_labeled(word)
```

Standardize a word into `letterpair` biletters ``(letter, occurrence number)``.

<a id="schubmult.combinatorics.indexed_forests.word_to_pairinj_labeled"></a>

#### word\_to\_pairinj\_labeled

```python
def word_to_pairinj_labeled(word)
```

Standardize a word into `ParallelInjLetter` biletters ``(letter, occurrence number)``.

<a id="schubmult.combinatorics.indexed_forests.omega_insertion"></a>

#### omega\_insertion

```python
def omega_insertion(
        word_of_pairs: tuple[letterpair,
                             ...]) -> tuple[LBS, DecLabeling] | None
```

Forest analogue of RSK: insert a word of biletters letter by letter, returning the insertion
forest ``P`` (an `LBS`) and the recording forest ``Q`` (a `DecLabeling`), or ``None`` if the
word is not insertable. Each new letter becomes a root that absorbs the neighboring trees
at ``primary - 1`` / ``primary + 1`` as children.

<a id="schubmult.combinatorics.indexed_forests.omega_reduced_word_from_labelings"></a>

#### omega\_reduced\_word\_from\_labelings

```python
def omega_reduced_word_from_labelings(P: LBS,
                                      Q: DecLabeling) -> tuple[int, ...]
```

Read the reduced-word order from an omega insertion pair (P, Q).

The decreasing labeling Q orders the surviving insertion events; sorting the
forest nodes by Q-label and reading primary labels from P gives the reduced
word represented by the omega insertion output.

<a id="schubmult.combinatorics.indexed_forests.omega_set_valued_compatibility"></a>

#### omega\_set\_valued\_compatibility

```python
def omega_set_valued_compatibility(word_of_pairs: tuple[letterpair, ...])
```

Check set-valued compatibility of a longer sequence via omega insertion.

  Returns a dict with:
      - ok: whether the sequence is set-valued compatible in WCGraph sense,
- omega_reduced_word: reduced word extracted from (P, Q),
- wc_reduced_word: reduced word from WCGraph reduction,
- set_sequence: reduced compatible set sequence from WCGraph.

<a id="schubmult.combinatorics.indexed_forests.omega_is_set_valued_compatible"></a>

#### omega\_is\_set\_valued\_compatible

```python
def omega_is_set_valued_compatible(
        word_of_pairs: tuple[letterpair, ...]) -> bool
```

Boolean convenience wrapper for omega_set_valued_compatibility.

<a id="schubmult.combinatorics.indexed_forests.omega_setvalued_insertion"></a>

#### omega\_setvalued\_insertion

```python
def omega_setvalued_insertion(word, compatible_sequence)
```

Omega insertion for potentially unreduced words with set-valued labels.

This mirrors the reduction logic used by WCGraph: scan a (word, compatible
sequence) pair left-to-right. If appending a letter would be non-reduced,
do not insert a new node; instead merge the compatible label into the set
attached to the corresponding reduced root. Otherwise perform ordinary
insertion on that letterpair.

Returns
-------
dict
    {
      "P": LBS,
      "Q": DecLabeling,
      "reduced_word": tuple[int, ...],
      "set_sequence": tuple[tuple[int, ...], ...],
      "node_sets": dict[int, tuple[int, ...]],
    }
where ``node_sets`` maps forest node index -> merged compatible-label set.

<a id="schubmult.combinatorics.indexed_forests.omega_setvalued_insertion_from_wcgraph"></a>

#### omega\_setvalued\_insertion\_from\_wcgraph

```python
def omega_setvalued_insertion_from_wcgraph(wc_graph)
```

Convenience wrapper: run set-valued omega insertion directly from WCGraph.

<a id="schubmult.combinatorics.indexed_forests.omega_park"></a>

#### omega\_park

```python
def omega_park(word)
```

Parking procedure Omega.

Cars 1, 2, ... arrive successively; car i prefers spot ``word[i]``.
If the preferred spot is empty, the car parks there. Otherwise the
preferred spot lies in a maximal interval [a, b] of occupied spots.
Let j < i be maximal with ``word[j] in [a, b]``; if ``word[i] >= word[j]``
the car parks at b+1, else at a-1.

Returns the set of occupied spots after all cars have parked.

<a id="schubmult.combinatorics.inversions_tableau"></a>

# schubmult.combinatorics.inversions\_tableau

Inversions tableaux: a root -> label-set assignment on the positive roots ``(i, j)``
(``i < j``) of a permutation, generalizing RC/WC graphs to a set-valued labeling.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau"></a>

## InversionsTableau Objects

```python
class InversionsTableau()
```

A dict-like assignment of label sets to roots ``(i, j)`` of a permutation, satisfying the
three axioms checked by `is_valid`. Build via the constructor (dict of root -> int-or-set),
`from_rc_graph`, or `from_wc_graph`; convert back via `to_rc_graph`/`to_wc_graph`.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.__init__"></a>

#### \_\_init\_\_

```python
def __init__(_dict, *_, **__)
```

Build from a mapping ``{(i, j): label_or_label_set}`` (roots to integer or set/frozenset labels).

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.iter_keys"></a>

#### iter\_keys

```python
def iter_keys(reverse=False)
```

Iterates in word order

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.iter_items"></a>

#### iter\_items

```python
def iter_items(reverse=False)
```

Iterates in word order

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc)
```

Build from an `RCGraph`: each left-to-right inversion root gets the row it's crossed at as its label.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.from_wc_graph"></a>

#### from\_wc\_graph

```python
@classmethod
def from_wc_graph(cls, wc)
```

Build from a `WCGraph`: replays its perm word/compatible sequence, tracking transported
roots and accumulating labels at each simple-root position.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.perm_word"></a>

#### perm\_word

```python
@cached_property
def perm_word()
```

A (not necessarily reduced) word recovering the tableau's permutation, built by repeatedly
peeling the largest label off a simple root and transporting the remaining roots.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.is_valid"></a>

#### is\_valid

```python
@property
def is_valid()
```

Whether the root/label assignment satisfies the three defining axioms:
(1) simple-root labels are bounded by the smaller index, (2) roots sharing a second
index have disjoint label sets, and (3) the label sets of ``(i,j)``, ``(j,k)``, ``(i,k)``
interleave consistently whenever all three roots are present.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.set_leq"></a>

#### set\_leq

```python
@staticmethod
def set_leq(set1, set2)
```

Check if set1 <= set2, i.e. max(set1) <= min(set2).

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.set_lt"></a>

#### set\_lt

```python
@staticmethod
def set_lt(set1, set2)
```

Check if set1 < set2, i.e. max(set1) < min(set2).

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.all_set_valued_inversions_tableaux"></a>

#### all\_set\_valued\_inversions\_tableaux

```python
@classmethod
def all_set_valued_inversions_tableaux(cls, perm, max_value=None)
```

Enumerate all set-valued inversions tableaux for ``perm``.

This enumerates candidate root-label assignments directly from the
defining axioms in :meth:`is_valid` (no WCGraph construction), then
filters by validity and target permutation.

The optional ``max_value`` bounds all labels by ``{1, ..., max_value}``.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.perm"></a>

#### perm

```python
@cached_property
def perm()
```

The permutation recovered from `perm_word` via the 0-Hecke (Demazure) product.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.compatible_sequence"></a>

#### compatible\_sequence

```python
@cached_property
def compatible_sequence()
```

The compatible sequence paired with `perm_word`: labels in weakly increasing order,
each repeated by the size of its label set.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.is_reduced"></a>

#### is\_reduced

```python
@cached_property
def is_reduced()
```

Whether ``perm_word`` is a reduced word for ``perm``.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.reduced_word"></a>

#### reduced\_word

```python
@cached_property
def reduced_word()
```

A reduced word for ``perm``, obtained from ``perm_word`` by dropping non-ascending steps
and folding their compatible-sequence entries into the surviving root's label set.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.is_set_valued"></a>

#### is\_set\_valued

```python
@property
def is_set_valued()
```

Whether any root has more than a single label (a genuinely set-valued tableau).

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph(length=None)
```

Convert a reduced (non-set-valued) tableau to an `RCGraph` via ``RCGraph.from_reduced_compatible``.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.to_wc_graph"></a>

#### to\_wc\_graph

```python
def to_wc_graph(length=None)
```

Convert to a `WCGraph` built directly from the root/label-set dictionary.

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x, y=None, *, beta=None, prop_beta=False)
```

Monomial contribution of this tableau to the (beta-deformed) Grothendieck polynomial:
``prod_v x[v] ** |labels at v|``, times a power of ``beta`` accounting for non-reduced excess.

<a id="schubmult.combinatorics.mbpd"></a>

# schubmult.combinatorics.mbpd

Implementation of the MBPD <-> RCP <-> WCGraph bijection from
``writing/mbpd.solve.tex`` (Theorem "T: main").

Geometry (paper convention):
  * n x n grid, rows top->bottom (1..n), columns left->right (1..n).
  * Pipes ENTER from the RIGHT border and EXIT to the SOUTH border.
  * No pipe enters from the top, no pipe leaves to the left.

Seven tiles, described by the set of directions in which the pipe connects
(N=up, E=right, S=down, W=left) together with a "marked" bit:

    B  blank        {}            heavy
    H  horizontal   {E, W}
    V  vertical     {N, S}
    P  plus/cross   {N, E, S, W}
    R  R-elbow      {E, S}
    J  J-elbow      {N, W}
    M  marked J     {N, W}        heavy

A tile is *heavy* iff it is B or M.

The bijection ``Phi = MBPD.phi`` and its inverse ``Psi = RCP.psi`` preserve the
associated permutation and weight.  ``WCGraph.to_mbpd`` / ``WCGraph.from_mbpd``
expose the composition with the trivial ``RCP <-> WCGraph`` repackaging.

<a id="schubmult.combinatorics.mbpd.tile_name"></a>

#### tile\_name

```python
def tile_name(conn: frozenset, marked: bool) -> str
```

Return the canonical tile name for a connection set / mark.

<a id="schubmult.combinatorics.mbpd.is_heavy"></a>

#### is\_heavy

```python
def is_heavy(conn: frozenset, marked: bool) -> bool
```

A tile is heavy iff it is blank or a marked J.

<a id="schubmult.combinatorics.mbpd.MBPD"></a>

## MBPD Objects

```python
class MBPD()
```

A marked bumpless pipedream on an ``n x n`` grid.

Internally we store, for every cell ``(i, j)`` (1-indexed), the connection
set ``conn[i][j]`` (a frozenset of ``N/E/S/W``) and a boolean ``marked``.

<a id="schubmult.combinatorics.mbpd.MBPD.__init__"></a>

#### \_\_init\_\_

```python
def __init__(n: int, conn, marked)
```

Store the ``n x n`` connection sets and marks as tuples (see `from_tiles` for the string form).

<a id="schubmult.combinatorics.mbpd.MBPD.conn"></a>

#### conn

```python
def conn(i: int, j: int) -> frozenset
```

Connection set of cell ``(i, j)`` (1-indexed).

<a id="schubmult.combinatorics.mbpd.MBPD.marked"></a>

#### marked

```python
def marked(i: int, j: int) -> bool
```

Whether cell ``(i, j)`` is marked.

<a id="schubmult.combinatorics.mbpd.MBPD.tile"></a>

#### tile

```python
def tile(i: int, j: int) -> str
```

Tile name (``B/H/V/P/R/J/M``) of cell ``(i, j)``.

<a id="schubmult.combinatorics.mbpd.MBPD.heavy"></a>

#### heavy

```python
def heavy(i: int, j: int) -> bool
```

Whether cell ``(i, j)`` is heavy.

<a id="schubmult.combinatorics.mbpd.MBPD.connects"></a>

#### connects

```python
def connects(i: int, j: int, d: str) -> bool
```

Whether the pipe in cell ``(i, j)`` connects in direction ``d``.

<a id="schubmult.combinatorics.mbpd.MBPD.from_tiles"></a>

#### from\_tiles

```python
@classmethod
def from_tiles(cls, grid) -> MBPD
```

Build from an ``n x n`` grid of tile-name strings.

<a id="schubmult.combinatorics.mbpd.MBPD.to_tiles"></a>

#### to\_tiles

```python
def to_tiles()
```

The ``n x n`` grid of tile-name strings.

<a id="schubmult.combinatorics.mbpd.MBPD.with_tile"></a>

#### with\_tile

```python
def with_tile(i: int, j: int, name: str) -> MBPD
```

Return a copy with cell ``(i, j)`` set to tile ``name``.

<a id="schubmult.combinatorics.mbpd.MBPD.rothe"></a>

#### rothe

```python
@classmethod
def rothe(cls, perm, n=None) -> MBPD
```

The Rothe MBPD ``D_w``: pipe ``i`` turns only at ``(i, w(i))``.

Connections (derived and verified against the paper's conventions):
  E : j >= w(i)                 (horizontal present / border on right)
  W : j > w(i)
  S : w^{-1}(j) <= i            (vertical present going down)
  N : w^{-1}(j) < i

<a id="schubmult.combinatorics.mbpd.MBPD.validity_errors"></a>

#### validity\_errors

```python
def validity_errors()
```

Return a list of human-readable validity problems (empty if valid).

<a id="schubmult.combinatorics.mbpd.MBPD.is_valid"></a>

#### is\_valid

```python
def is_valid() -> bool
```

Whether `validity_errors` is empty.

<a id="schubmult.combinatorics.mbpd.MBPD.perm"></a>

#### perm

```python
def perm()
```

Associated permutation ``w`` where the pipe entering the right of
row ``i`` leaves the bottom of column ``w(i)``.

Double crossings between an already-crossed pair are ignored (the
Demazure convention): the *first* time a pair of pipes meets at a ``P``
tile they cross, every later meeting is a *bump* (the ``P`` acts as
superimposed elbows ``R``+``J``, i.e. ``E<->S`` and ``N<->W``).

We compute the routing by a fixed point.  Start assuming every ``P`` is
a genuine crossing (straight through) and trace all pipes.  Whenever a
pair of pipes meets at more than one ``P`` tile we mark all but their
first meeting as bumps, then re-trace.  This is repeated until the set
of bump tiles stabilises; the exit columns then give ``w``.

<a id="schubmult.combinatorics.mbpd.MBPD.weight"></a>

#### weight

```python
def weight()
```

``wt(D) = (m_1, ..., m_n)`` with ``m_i`` = [`heavy`](#schubmult.combinatorics.mbpd.MBPD.heavy) tiles in row ``i``.

<a id="schubmult.combinatorics.mbpd.MBPD.num_heavy"></a>

#### num\_heavy

```python
def num_heavy() -> int
```

Total number of heavy tiles.

<a id="schubmult.combinatorics.mbpd.MBPD.heavy_cells"></a>

#### heavy\_cells

```python
def heavy_cells()
```

Positions ``(i, j)`` of all heavy tiles in row-major order.

<a id="schubmult.combinatorics.mbpd.MBPD.is_pipe_segment"></a>

#### is\_pipe\_segment

```python
def is_pipe_segment(r: int, b: int, c: int) -> bool
```

``D_{r,[b,c]}`` is a pipe segment: a connected horizontal run.

For ``b == c`` a single non-blank tile qualifies.  For ``b < c`` the
interior ``D_{r,[b+1,c-1]}`` must be all ``H``/``P`` tiles and the run
must actually connect: ``(r,b)`` connects east and ``(r,c)`` connects
west (so the paper's inference that the endpoints connect toward the
interior holds).

<a id="schubmult.combinatorics.mbpd.MBPD.rj_subsequence"></a>

#### rj\_subsequence

```python
def rj_subsequence(r: int, lo: int, hi: int)
```

The RJ subsequence (list of 'R'/'J') of light tiles in ``[lo,hi]``.

Returns ``None`` if any tile in the range is heavy (not a light seq).

<a id="schubmult.combinatorics.mbpd.MBPD.light_seq_type"></a>

#### light\_seq\_type

```python
def light_seq_type(r: int, lo: int, hi: int)
```

Classify the light sequence ``D_{r,[lo,hi]}`` as one of
``'paired'``, ``'J'``, ``'R'``, ``'JR'`` (or ``None`` if it is not a
light sequence).  Empty range counts as ``'paired'``.

<a id="schubmult.combinatorics.mbpd.MBPD.is_doublecross"></a>

#### is\_doublecross

```python
def is_doublecross(r: int, b: int, d: int) -> bool
```

``D_{[r,r+1],[b,d]}`` is a doublecross: both rows are pipe segments,
``D_{r,b}=R`` and ``D_{r+1,d}=J``.

<a id="schubmult.combinatorics.mbpd.MBPD.admits_droop"></a>

#### admits\_droop

```python
def admits_droop(r: int, b: int, d: int) -> bool
```

Whether the ``(r, [b, d])``-droop (moving a pipe segment from row ``r`` down to row ``r+1`` over
columns ``b..d``) is admitted by the local tile configuration.

<a id="schubmult.combinatorics.mbpd.MBPD.admits_undroop"></a>

#### admits\_undroop

```python
def admits_undroop(r: int, b: int, d: int) -> bool
```

Whether the inverse of the ``(r, [b, d])``-droop is admitted.

<a id="schubmult.combinatorics.mbpd.MBPD.droop"></a>

#### droop

```python
def droop(r: int, b: int, d: int) -> MBPD
```

Apply the ``(r, [b, d])``-droop (raises if not admitted).

<a id="schubmult.combinatorics.mbpd.MBPD.undroop"></a>

#### undroop

```python
def undroop(r: int, b: int, d: int) -> MBPD
```

Apply the ``(r, [b, d])``-undroop (raises if not admitted).

<a id="schubmult.combinatorics.mbpd.MBPD.is_f_target"></a>

#### is\_f\_target

```python
def is_f_target(r: int, c: int) -> bool
```

Whether the heavy tile at ``(r, c)`` is a target of the ``f`` move of the bijection ``Phi``.

<a id="schubmult.combinatorics.mbpd.MBPD.is_fstar_target"></a>

#### is\_fstar\_target

```python
def is_fstar_target(r: int, c: int) -> bool
```

Whether the heavy tile at ``(r, c)`` is a target of the ``f*`` (terminal) move of ``Phi``.

<a id="schubmult.combinatorics.mbpd.MBPD.is_F_target"></a>

#### is\_F\_target

```python
def is_F_target(r: int, c: int) -> bool
```

Whether ``(r, c)`` is an ``f`` or ``f*`` target.

<a id="schubmult.combinatorics.mbpd.MBPD.max_F_target"></a>

#### max\_F\_target

```python
def max_F_target()
```

Bottommost then rightmost heavy tile, or ``None`` for ``D_id``.

<a id="schubmult.combinatorics.mbpd.MBPD.is_F_terminal"></a>

#### is\_F\_terminal

```python
def is_F_terminal() -> bool
```

Whether the next ``Phi`` step is terminal (no heavy tiles, or the max target is an ``f*`` target).

<a id="schubmult.combinatorics.mbpd.MBPD.F_target_info"></a>

#### F\_target\_info

```python
def F_target_info(r: int, c: int) -> dict
```

Describe the ``Phi`` move at the target ``(r, c)``: its kind (``f``/``fstar``), the resulting
MBPD, and the biletter emitted.

<a id="schubmult.combinatorics.mbpd.MBPD.f_move"></a>

#### f\_move

```python
def f_move(r: int, c: int) -> MBPD
```

Apply the ``F``-move at the ``F``-target ``(r,c)`` (paper
\S "do F move").

<a id="schubmult.combinatorics.mbpd.MBPD.e_target"></a>

#### e\_target

```python
def e_target(r: int)
```

Determine the ``E``-target ``(r+1,c)`` in the two-row window with top
row ``r`` (``1 <= r < n``).  Returns a dict with keys
``kind`` (``'e'`` or ``'estar'``), ``c``, ``cprime`` or ``None`` if no
``E``-target exists for this row.

<a id="schubmult.combinatorics.mbpd.MBPD.E_target_info"></a>

#### E\_target\_info

```python
def E_target_info(r: int) -> dict
```

Full case analysis for the ``E``-move at row ``r``: the right
trichotomy (Initial/Plus/NPlus) with ``rho``, and the left trichotomy
(Straight/DCross/Leftturn) with ``lambda``.

<a id="schubmult.combinatorics.mbpd.MBPD.e_move"></a>

#### e\_move

```python
def e_move(r: int) -> MBPD
```

Apply the ``E``-move ``E_r`` at row ``r`` (paper \S "the E moves"):
* if ``(r+1,c)`` is ``M``, unmark it;
* in Cases Straight and DCross do the ``(r,[lambda,rho])``-droop;
* if ``(r,lambda)`` is ``J``, mark it.

<a id="schubmult.combinatorics.mbpd.MBPD.phi"></a>

#### phi

```python
def phi() -> RCP
```

The bijection ``Phi`` of the paper (Theorem "T: main"), realized by
the row-pop recursion:

  * ``Phi(D_id) = ()``
  * F-terminal:    ``Phi(D) = (i,i) . Phi(f*_i(D))``
  * F-nonterminal: ``Phi(D) = down( Phi(f_i(D)) )``

where ``i`` is the row of the maximum ``F``-target and ``down`` replaces
the first biletter ``(i',a)`` by ``(i'-1,a)`` (RRCP2).

<a id="schubmult.combinatorics.mbpd.MBPD.to_bpd"></a>

#### to\_bpd

```python
def to_bpd()
```

Return ``(bpd, marks)``: the underlying :class:`BPD` and the frozenset
of ``(i, j)`` cells (1-indexed) that carry a mark (``M`` tiles).

<a id="schubmult.combinatorics.mbpd.MBPD.from_bpd"></a>

#### from\_bpd

```python
@classmethod
def from_bpd(cls, bpd, marks=frozenset()) -> MBPD
```

Rebuild an :class:`MBPD` from a :class:`BPD` and a set of marked
``(i, j)`` cells.  The BPD is squared up to ``max(rows, len(perm))`` first
(Rothe completion), and a mark is only reinstated where the tile is a
``J`` up-elbow.

<a id="schubmult.combinatorics.mbpd.MBPD.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
def zero_out_last_row(rows: int) -> MBPD
```

Zero out below ``rows`` by resizing the underlying BPD down to ``rows``
rows (Weigandt/Lascoux transition) and re-deducing the marks that survive.

Weight is preserved; the permutation changes according to the transition.

<a id="schubmult.combinatorics.mbpd.RCP"></a>

## RCP Objects

```python
@dataclass(frozen=True)
class RCP()
```

A reverse compatible pair: a tuple of biletters ``(i, a)`` with
``1 <= i <= a < n`` (``n`` is the ambient size), strictly *decreasing* in
the order ``(i1,a1) > (i2,a2)`` iff ``i1 > i2`` or (``i1 == i2`` and
``a1 < a2``).

<a id="schubmult.combinatorics.mbpd.RCP.biletters"></a>

#### biletters

tuple of (i, a)

<a id="schubmult.combinatorics.mbpd.RCP.is_valid"></a>

#### is\_valid

```python
def is_valid() -> bool
```

Whether every biletter satisfies ``1 <= i <= a < n`` and the sequence is strictly decreasing.

<a id="schubmult.combinatorics.mbpd.RCP.weight"></a>

#### weight

```python
def weight()
```

Number of biletters with each first coordinate ``i``, as a length-``n`` tuple.

<a id="schubmult.combinatorics.mbpd.RCP.perm"></a>

#### perm

```python
def perm()
```

``w(B) = s_{a_l} * ... * s_{a_1}`` (Demazure product, subscripts
decreasing = reversed biletter order).

<a id="schubmult.combinatorics.mbpd.RCP.to_wcgraph"></a>

#### to\_wcgraph

```python
def to_wcgraph()
```

Package the RCP as a :class:`WCGraph`.

``WCGraph.from_word_compatible`` wants a *weakly increasing* compatible
sequence.  Reversing the biletter list turns the RCP's decreasing order
into weakly increasing ``i`` with strictly decreasing ``a`` within a
block -- exactly the WCGraph compatibility condition.

<a id="schubmult.combinatorics.mbpd.RCP.from_wcgraph"></a>

#### from\_wcgraph

```python
@classmethod
def from_wcgraph(cls, wc, n=None)
```

Inverse of :meth:`to_wcgraph`.

<a id="schubmult.combinatorics.mbpd.RCP.to_pd_crossings"></a>

#### to\_pd\_crossings

```python
def to_pd_crossings()
```

Remark "R:CP and PD": biletter ``(i, a)`` -> crossing at
``(i, a - i + 1)`` in the ordinary pipedream picture.

<a id="schubmult.combinatorics.mbpd.RCP.psi"></a>

#### psi

```python
def psi() -> MBPD
```

The inverse bijection ``Psi = Phi^{-1}`` (Theorem "T: main"),
realized by row-unpop.

``Phi(D) = (i,a) . Phi(nabla_r D)`` where the maximum ``F``-target of
``D`` is in row ``i`` and ``nabla_r D = f*_a f_{a-1} ... f_i (D)`` is
obtained by climbing rows ``i, i+1, ..., a`` (the last, at row ``a``,
being the ``f*``-move).  Inverting: given the first biletter ``(i,a)``,
first rebuild ``N = nabla_r D = Psi(rest)``, then apply
``D = e_i(e_{i+1}( ... e_{a-1}( e*_a(N) ) ... ))``.

<a id="schubmult.combinatorics.nilplactic"></a>

# schubmult.combinatorics.nilplactic

`NilPlactic`: the nilCoxeter/nilHecke analogue of `Plactic`, used for Edelman-Greene tableaux
(skew shapes over increasing/decreasing-adjacent tableau rules) and their reduced-word data.

<a id="schubmult.combinatorics.nilplactic.NilPlactic"></a>

## NilPlactic Objects

```python
class NilPlactic(Plactic)
```

A nilCoxeter/nilHecke skew tableau: like `Plactic` but under the Edelman-Greene insertion
rule (increasing rows and columns, no repeated entries).

<a id="schubmult.combinatorics.nilplactic.NilPlactic.all_skew_ed_tableaux"></a>

#### all\_skew\_ed\_tableaux

```python
@classmethod
@cache
def all_skew_ed_tableaux(cls, outer_shape, bruhat_perm, inner_shape=None)
```

Return all skew Edelman–Greene tableaux of skew shape outer_shape/inner_shape
(represented as row lengths) that are increasing (strictly increasing rows
left-to-right and columns top-to-bottom) and whose row-reading word (E-G
row word: read rows bottom-to-top, left-to-right, omitting zeros) has
Permutation.ref_product(...) less than or equal to ``bruhat_perm`` in
Bruhat order.

Representation details:

- outer_shape: sequence of row lengths (one int per row).
- inner_shape: optional sequence of left offsets per row (defaults to zeros).
- The returned tableaux are NilPlactic instances whose rows are left-aligned
  of length ``outer_shape[r]``, with the first ``inner_shape[r]`` entries set
  to 0 to represent the inner (removed) cells of the skew shape.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(row, col)
```

Perform an upward K-theoretic jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new Plactic tableau.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(row, col)
```

Perform a K-theoretic jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new NilPlactic tableau.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.ed_insert"></a>

#### ed\_insert

```python
def ed_insert(*letters)
```

Insert a letter/entry into this NilPlactic tableau and return a new Plactic.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.ed_insert_rsk"></a>

#### ed\_insert\_rsk

```python
@classmethod
def ed_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this NilPlactic tableau and return a new Plactic.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.ed_column_insert_rsk"></a>

#### ed\_column\_insert\_rsk

```python
@classmethod
def ed_column_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this NilPlactic tableau and return a new Plactic.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other)
```

Plactic product: insert entries of `other` in row-reading order
(top-to-bottom, left-to-right) into a copy of self.

<a id="schubmult.combinatorics.permutation"></a>

# schubmult.combinatorics.permutation

The core indexing object used throughout ``schubmult``: finite permutations.

A `Permutation` wraps a 1-indexed array (``self[i]`` is the value at 0-indexed
position ``i``, extended by the identity ``self[i] = i + 1`` past its window)
and is immutable/hashable/cached (equal permutations of different window
lengths, e.g. ``[2, 1]`` and ``[2, 1, 3]``, compare and hash equal).
Multiplication ``p * q`` is composition of the underlying functions,
``(p * q)[i] = p[q[i] - 1]``; ``~p`` is the inverse. Most of the combinatorics
library is built from a `Permutation`'s Lehmer code (``code``/``trimcode``),
inversions, and Bruhat/weak order.

<a id="schubmult.combinatorics.permutation.Permutation"></a>

## Permutation Objects

```python
class Permutation(Printable)
```

A finite permutation, stored as a 1-indexed array and extended by the identity.

Construct from array form, e.g. ``Permutation([2, 1, 3])``, or via
``uncode(lehmer_code)``. Supports composition (``*``), inversion (``~``),
Bruhat order (``<=``/``bruhat_leq``), and iteration over the window
(``list(perm)``, ``perm[i]`` 0-indexed).

<a id="schubmult.combinatorics.permutation.Permutation.apply"></a>

#### apply

```python
def apply(arr)
```

Permute the elements of ``arr`` by this permutation: return ``[arr[self[i] - 1] for i in range(len(arr))]``.

<a id="schubmult.combinatorics.permutation.Permutation.is_reducible"></a>

#### is\_reducible

```python
@property
def is_reducible()
```

Whether ``self`` splits as a direct sum of two smaller permutations across some fixed point.

<a id="schubmult.combinatorics.permutation.Permutation.reduce"></a>

#### reduce

```python
def reduce(start_spot=1, strict=False)
```

Split ``self`` at a fixed point into ``(perm1, perm2)`` on disjoint blocks of values, or ``None`` if not reducible.

<a id="schubmult.combinatorics.permutation.Permutation.act_root"></a>

#### act\_root

```python
def act_root(a, b)
```

Image of the root/pair ``(a, b)`` (1-indexed positions) under ``self``: ``(self[a-1], self[b-1])``.

<a id="schubmult.combinatorics.permutation.Permutation.one_dominates"></a>

#### one\_dominates

```python
def one_dominates(other)
```

Whether ``self`` one-step dominates ``other`` (one level of the recursive ``dominates`` test).

<a id="schubmult.combinatorics.permutation.Permutation.dominates"></a>

#### dominates

```python
def dominates(other)
```

Whether ``self`` dominates ``other`` in the sense used by the dual Pieri / positivity rules.

<a id="schubmult.combinatorics.permutation.Permutation.antiperm"></a>

#### antiperm

```python
@property
def antiperm()
```

Conjugate of ``self`` by the longest element ``w0``: ``w0 * self * w0``.

<a id="schubmult.combinatorics.permutation.Permutation.weight_coset_decomp"></a>

#### weight\_coset\_decomp

```python
def weight_coset_decomp(dominant_weight)
```

Coset decomposition of ``self`` with respect to the stabilizer of ``dominant_weight``; see ``coset_decomp``.

<a id="schubmult.combinatorics.permutation.Permutation.min_of_weight_coset"></a>

#### min\_of\_weight\_coset

```python
def min_of_weight_coset(dominant_weight)
```

Minimal-length coset representative of ``self`` in the stabilizer of ``dominant_weight``.

<a id="schubmult.combinatorics.permutation.Permutation.max_of_weight_coset"></a>

#### max\_of\_weight\_coset

```python
def max_of_weight_coset(dominant_weight)
```

Maximal-length coset representative of ``self`` in the stabilizer of ``dominant_weight``.

<a id="schubmult.combinatorics.permutation.Permutation.does_demazure_crystal_tensor_decompose"></a>

#### does\_demazure\_crystal\_tensor\_decompose

```python
@staticmethod
def does_demazure_crystal_tensor_decompose(dominant_weight1, low_perm1,
                                           dominant_weight2, low_perm2)
```

Whether the tensor product of the two Demazure crystals (given by their dominant weights and
lowest-weight sorting permutations) decomposes as their `does_demazure_crystal_tensor_decompose`
criterion predicts: the minimal coset representative of ``low_perm1`` in weight ``dominant_weight1``
lies below the longest element of the left descents of the maximal coset representative of
``low_perm2`` in weight ``dominant_weight2``, in Bruhat order.

<a id="schubmult.combinatorics.permutation.Permutation.coset_decomp"></a>

#### coset\_decomp

```python
def coset_decomp(*descs)
```

Decompose ``self = reduced_perm * w_J`` where ``w_J`` lies in the parabolic subgroup generated
by the simple reflections at 1-indexed positions ``descs`` and ``reduced_perm`` is the minimal-length
coset representative (no descents in ``descs``).

**Returns**:

- `tuple` - ``(reduced_perm, w_J)``.

<a id="schubmult.combinatorics.permutation.Permutation.min_coset_rep"></a>

#### min\_coset\_rep

```python
def min_coset_rep(*descs)
```

Minimal-length coset representative of ``self`` for the parabolic subgroup generated at ``descs``.

<a id="schubmult.combinatorics.permutation.Permutation.max_coset_rep"></a>

#### max\_coset\_rep

```python
def max_coset_rep(*descs)
```

Maximal-length coset representative of ``self`` for the parabolic subgroup generated at ``descs``.

<a id="schubmult.combinatorics.permutation.Permutation.longest_element"></a>

#### longest\_element

```python
@classmethod
def longest_element(cls, *descs)
```

Longest element of the parabolic subgroup generated by the simple reflections at 1-indexed positions ``descs``.

<a id="schubmult.combinatorics.permutation.Permutation.w0"></a>

#### w0

```python
@classmethod
def w0(cls, n)
```

The longest element of the symmetric group on ``n`` letters: ``[n, n-1, ..., 1]``.

<a id="schubmult.combinatorics.permutation.Permutation.all_permutations"></a>

#### all\_permutations

```python
@classmethod
@cache
def all_permutations(cls, n)
```

All permutations of ``{1, ..., n}``, as a list of `Permutation`.

<a id="schubmult.combinatorics.permutation.Permutation.parabolic_reduce"></a>

#### parabolic\_reduce

```python
def parabolic_reduce(*descs)
```

Decompose ``self = reduced_perm * w_J`` where ``w_J`` is generated by the simple reflections
*not* in ``descs`` and ``reduced_perm`` has no descents outside ``descs``.

**Returns**:

- `tuple` - ``(reduced_perm, w_J)``.

<a id="schubmult.combinatorics.permutation.Permutation.__truediv__"></a>

#### \_\_truediv\_\_

```python
def __truediv__(other)
```

Returns a tuple of (self, other) if other is a permutation. Intended for skew elements.

<a id="schubmult.combinatorics.permutation.Permutation.fixers"></a>

#### fixers

```python
@staticmethod
def fixers(dominant_weight)
```

1-indexed positions ``i`` where ``dominant_weight[i-1] == dominant_weight[i]``
(the simple reflections fixing ``dominant_weight``, i.e. generating its stabilizer).

<a id="schubmult.combinatorics.permutation.Permutation.ref_product"></a>

#### ref\_product

```python
@classmethod
def ref_product(cls, *args)
```

Product of the simple reflections ``s_a`` for ``a`` in ``args``, applied left to right.

<a id="schubmult.combinatorics.permutation.Permutation.hecke_ref_product"></a>

#### hecke\_ref\_product

```python
@classmethod
def hecke_ref_product(cls, *args)
```

Like ``ref_product``, but each ``s_a`` is applied only if it is a Bruhat ascent
(the 0-Hecke / Demazure product of the simple reflections).

<a id="schubmult.combinatorics.permutation.Permutation.code_word"></a>

#### code\_word

```python
@property
def code_word()
```

A canonical reduced word for ``self``, read off from ``trimcode``.

<a id="schubmult.combinatorics.permutation.Permutation.inverse_code_word"></a>

#### inverse\_code\_word

```python
@property
def inverse_code_word()
```

A canonical reduced word for ``~self``, read off from ``(~self).trimcode``.

<a id="schubmult.combinatorics.permutation.Permutation.root_swap"></a>

#### root\_swap

```python
def root_swap(root)
```

Multiply ``self`` by the reflection swapping the 1-indexed positions ``root = (a, b)``.

<a id="schubmult.combinatorics.permutation.Permutation.reflection"></a>

#### reflection

```python
@classmethod
def reflection(cls, root)
```

The transposition swapping the 1-indexed positions ``root = (a, b)``.

<a id="schubmult.combinatorics.permutation.Permutation.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(index, word=None)
```

The positive root sent to a negative root by the ``index``-th letter of ``word``
(default ``self.code_word``), read from the right (post-multiplied by the remaining suffix).

<a id="schubmult.combinatorics.permutation.Permutation.left_root_at"></a>

#### left\_root\_at

```python
def left_root_at(index, word=None)
```

Like ``right_root_at``, but the root is transported by the prefix of ``word`` before ``index``.

<a id="schubmult.combinatorics.permutation.Permutation.all_reduced_words"></a>

#### all\_reduced\_words

```python
@cache
def all_reduced_words()
```

All reduced words of `self`, by peeling descents.

<a id="schubmult.combinatorics.permutation.Permutation.all_reduced_subwords"></a>

#### all\_reduced\_subwords

```python
@staticmethod
def all_reduced_subwords(word)
```

All reduced subwords of `self`, by peeling descents.

<a id="schubmult.combinatorics.permutation.Permutation.all_subwords"></a>

#### all\_subwords

```python
@staticmethod
def all_subwords(word)
```

All subwords of `self`, by peeling descents.

<a id="schubmult.combinatorics.permutation.Permutation.commutation_class_of"></a>

#### commutation\_class\_of

```python
@staticmethod
def commutation_class_of(word)
```

All words obtainable from ``word`` by commuting adjacent far-apart letters (``|a - b| >= 2``).

<a id="schubmult.combinatorics.permutation.Permutation.forest_class_of"></a>

#### forest\_class\_of

```python
@staticmethod
def forest_class_of(word)
```

Words reachable from ``word`` by commutation moves that also preserve the indexed-forest
insertion structure (``omega_insertion``); a refinement of ``commutation_class_of``.

<a id="schubmult.combinatorics.permutation.Permutation.code_index_of_index"></a>

#### code\_index\_of\_index

```python
def code_index_of_index(index)
```

The position in ``trimcode`` whose block of ``code_word`` letters contains position ``index``.

<a id="schubmult.combinatorics.permutation.Permutation.cycle"></a>

#### cycle

```python
@staticmethod
def cycle(p, q)
```

Construct the cycle permutation used elsewhere in the code.
Kept as a staticmethod on Permutation for call sites like Permutation.cycle(p,q).

<a id="schubmult.combinatorics.permutation.Permutation.sorting_perm"></a>

#### sorting\_perm

```python
@classmethod
def sorting_perm(cls, itera, reverse=False)
```

The permutation that sorts ``itera`` into (by default) increasing order.

<a id="schubmult.combinatorics.permutation.Permutation.right_act"></a>

#### right\_act

```python
def right_act(lst)
```

Permute the entries of ``lst`` by this permutation, preserving ``lst``'s type (list stays a list).

<a id="schubmult.combinatorics.permutation.Permutation.bruhat_leq"></a>

#### bruhat\_leq

```python
def bruhat_leq(perm, perm2)
```

Whether ``perm <= perm2`` in Bruhat order (equivalently, ``perm``'s tableau criterion
against ``perm2`` on every prefix of their windows).

<a id="schubmult.combinatorics.permutation.Permutation.from_code"></a>

#### from\_code

```python
@classmethod
def from_code(cls, cd)
```

Alias for ``uncode(cd)``: the permutation with Lehmer code ``cd``.

<a id="schubmult.combinatorics.permutation.Permutation.has_pattern"></a>

#### has\_pattern

```python
def has_pattern(pattern)
```

Whether ``self`` contains ``pattern`` (a plain list, not necessarily reduced) as a pattern:
some subsequence of ``self``'s window order-isomorphic to ``pattern``.

<a id="schubmult.combinatorics.permutation.Permutation.zero_indexed_descents"></a>

#### zero\_indexed\_descents

```python
def zero_indexed_descents()
```

0-indexed descent positions: ``i`` such that ``self[i] > self[i+1]``.

<a id="schubmult.combinatorics.permutation.Permutation.descents"></a>

#### descents

```python
def descents(zero_indexed=True)
```

Descent positions of ``self``, 0-indexed by default or 1-indexed if ``zero_indexed=False``.

<a id="schubmult.combinatorics.permutation.Permutation.get_cycles"></a>

#### get\_cycles

```python
def get_cycles(sort_min=False)
```

Cycle decomposition of ``self`` as a list of tuples; ``sort_min`` rotates each cycle to start
at its minimum element instead of the sympy convention.

<a id="schubmult.combinatorics.permutation.Permutation.from_cycles"></a>

#### from\_cycles

```python
@classmethod
def from_cycles(cls, cycle_iter)
```

Build a permutation from a cycle decomposition (sequence of cycles, each a sequence of 1-indexed values).

<a id="schubmult.combinatorics.permutation.Permutation.code"></a>

#### code

```python
@property
def code()
```

Lehmer code of ``self``: ``code[i]`` counts ``j > i`` with ``self[i] > self[j]``.

<a id="schubmult.combinatorics.permutation.Permutation.graph"></a>

#### graph

```python
@property
def graph()
```

The permutation matrix support as a set of 1-indexed pairs ``{(i, self[i])}``.

<a id="schubmult.combinatorics.permutation.Permutation.reduced_with"></a>

#### reduced\_with

```python
def reduced_with(other)
```

Whether ``self * other`` is length-additive: ``inv(self) + inv(other) == inv(self * other)``.

<a id="schubmult.combinatorics.permutation.Permutation.shape"></a>

#### shape

```python
@property
def shape()
```

The Lehmer code sorted into weakly decreasing order (a partition).

<a id="schubmult.combinatorics.permutation.Permutation.inversion_set"></a>

#### inversion\_set

```python
@cached_property
def inversion_set()
```

Set of inversions ``(i, j)`` (1-indexed, ``i < j``) with ``self[i-1] > self[j-1]``.

<a id="schubmult.combinatorics.permutation.Permutation.weak_order_leq"></a>

#### weak\_order\_leq

```python
def weak_order_leq(other)
```

Whether ``self <= other`` in left weak order (``self``'s inversion set is a subset of ``other``'s).

<a id="schubmult.combinatorics.permutation.Permutation.weak_order_meet"></a>

#### weak\_order\_meet

```python
def weak_order_meet(other)
```

Meet (greatest lower bound) of ``self`` and ``other`` in left weak order.

<a id="schubmult.combinatorics.permutation.Permutation.weak_order_join"></a>

#### weak\_order\_join

```python
def weak_order_join(other)
```

Join (least upper bound) of ``self`` and ``other`` in left weak order, via ``w0``-duality with the meet.

<a id="schubmult.combinatorics.permutation.Permutation.diagram"></a>

#### diagram

```python
@property
def diagram()
```

The Rothe diagram of ``self``: cells ``(i, j)`` with ``self[i-1] > j`` and ``(~self)[j-1] > i``.

<a id="schubmult.combinatorics.permutation.Permutation.rothe_diagram"></a>

#### rothe\_diagram

```python
@property
def rothe_diagram()
```

The graph of ``self`` as a set of 1-indexed pairs (see ``diagram`` for the Rothe diagram cells).

<a id="schubmult.combinatorics.permutation.Permutation.max_descent"></a>

#### max\_descent

```python
@cached_property
def max_descent()
```

Number of entries in ``trimcode`` (one past the last nonzero code entry).

<a id="schubmult.combinatorics.permutation.Permutation.maximal_corner"></a>

#### maximal\_corner

```python
@property
def maximal_corner()
```

The maximal corner ``(maxd, end_spot)`` of ``self``'s diagram, used by ``pivots``/``pivot_transition``.

<a id="schubmult.combinatorics.permutation.Permutation.from_partial"></a>

#### from\_partial

```python
@classmethod
def from_partial(cls, partial_perm)
```

Complete a partial assignment (a list with some ``None`` entries) to a full permutation,
filling the gaps with the missing values in increasing order.

<a id="schubmult.combinatorics.permutation.Permutation.pivots"></a>

#### pivots

```python
@cache
def pivots(a=None, b=None)
```

Return the set of pivot positions for a maximal corner (a,b).

<a id="schubmult.combinatorics.permutation.Permutation.pivot_transition"></a>

#### pivot\_transition

```python
def pivot_transition(pivot_set)
```

Grothendieck transition for a given pivot set at the maximal corner. Returns the resulting permutation.

<a id="schubmult.combinatorics.permutation.Permutation.pad_code"></a>

#### pad\_code

```python
def pad_code(length)
```

``trimcode`` padded with trailing zeros to ``length``.

<a id="schubmult.combinatorics.permutation.Permutation.trimcode"></a>

#### trimcode

```python
@cached_property
def trimcode()
```

Lehmer code truncated to drop trailing zeros (length equals the last descent position).

<a id="schubmult.combinatorics.permutation.Permutation.mul_dominant"></a>

#### mul\_dominant

```python
def mul_dominant()
```

Left factor of ``self`` in its decomposition against the minimal dominant permutation above it.

<a id="schubmult.combinatorics.permutation.Permutation.strict_mul_dominant"></a>

#### strict\_mul\_dominant

```python
def strict_mul_dominant(size=None)
```

Variant of ``mul_dominant`` built from the strict theta (strictly decreasing dominant code).

<a id="schubmult.combinatorics.permutation.Permutation.shiftup"></a>

#### shiftup

```python
def shiftup(k)
```

``self`` with its Lehmer code shifted right by ``k`` (``k`` leading zero code entries prepended).

<a id="schubmult.combinatorics.permutation.Permutation.inv"></a>

#### inv

```python
@cached_property
def inv()
```

Length of ``self``: the number of inversions, i.e. ``sum(self.code)``.

<a id="schubmult.combinatorics.permutation.Permutation.is_dominant"></a>

#### is\_dominant

```python
@property
def is_dominant()
```

Whether ``self`` equals the minimal dominant permutation above it (its code is already weakly decreasing).

<a id="schubmult.combinatorics.permutation.Permutation.is_strict_dominant"></a>

#### is\_strict\_dominant

```python
@property
def is_strict_dominant()
```

Whether ``self`` is dominant with a strictly decreasing ``trimcode``.

<a id="schubmult.combinatorics.permutation.Permutation.is_vexillary"></a>

#### is\_vexillary

```python
@property
def is_vexillary()
```

Whether ``self`` avoids the pattern ``2143`` (equivalently, its Schubert polynomial is a
single Schur polynomial in the ``trimcode``-shape).

<a id="schubmult.combinatorics.permutation.Permutation.swap"></a>

#### swap

```python
def swap(i, j)
```

Multiply ``self`` by the transposition of 0-indexed positions ``i`` and ``j`` (window extended as needed).

<a id="schubmult.combinatorics.permutation.Permutation.rslice"></a>

#### rslice

```python
def rslice(start, stop)
```

Window values at 0-indexed positions ``start`` (inclusive) to ``stop`` (exclusive), extended by fixed points.

<a id="schubmult.combinatorics.permutation.Permutation.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

Demazure product

<a id="schubmult.combinatorics.permutation.Permutation.pattern_at"></a>

#### pattern\_at

```python
def pattern_at(*indices)
```

The (inverse sorting) pattern induced by ``self`` on the given 1-indexed ``indices``.

<a id="schubmult.combinatorics.permutation.Permutation.minimal_dominant_above"></a>

#### minimal\_dominant\_above

```python
def minimal_dominant_above()
```

The minimal dominant permutation ``>= self`` in Bruhat order: ``uncode(self.theta())``.

<a id="schubmult.combinatorics.permutation.Permutation.foundational_root"></a>

#### foundational\_root

```python
@property
def foundational_root()
```

The pivot pair ``(k, mx)`` marking ``self``'s last descent block and the last position dropping below it.

<a id="schubmult.combinatorics.permutation.Permutation.theta"></a>

#### theta

```python
def theta()
```

Dominant (weakly decreasing) sequence bounding ``self``'s code, used throughout the v-path
multiplication algorithms (see ``schubmult.mult``).

<a id="schubmult.combinatorics.permutation.Permutation.medium_theta"></a>

#### medium\_theta

```python
def medium_theta()
```

Variant of ``theta`` used by the "fast"/merged-layer multiplication kernels.

<a id="schubmult.combinatorics.permutation.Permutation.strict_theta"></a>

#### strict\_theta

```python
def strict_theta()
```

Variant of ``theta`` with strictly decreasing entries (no repeated nonzero layers).

<a id="schubmult.combinatorics.permutation.Permutation.maximal_sortable_below"></a>

#### maximal\_sortable\_below

```python
def maximal_sortable_below()
```

The maximal "sortable" permutation ``<= self`` (code has no gap of more than 1 between
consecutive entries), used by ``mul_sortable``.

<a id="schubmult.combinatorics.permutation.Permutation.mul_sortable"></a>

#### mul\_sortable

```python
def mul_sortable()
```

Left factor of ``self`` against its ``maximal_sortable_below`` decomposition (dual to ``mul_dominant``).

<a id="schubmult.combinatorics.permutation.uncode"></a>

#### uncode

```python
def uncode(cd)
```

The permutation whose Lehmer code is ``cd`` (a list of nonnegative integers).

<a id="schubmult.combinatorics.permutation.permtrim"></a>

#### permtrim

```python
def permtrim(perm)
```

Normalize ``perm`` (a plain array) into a `Permutation` (trims trailing fixed points).

<a id="schubmult.combinatorics.permutation.phi1"></a>

#### phi1

```python
def phi1(u)
```

Drop the first entry of ``(~u).code`` and re-invert: one step of the ``dominates`` recursion.

<a id="schubmult.combinatorics.permutation.split_perms"></a>

#### split\_perms

```python
def split_perms(perms)
```

Split each permutation in ``perms`` (after the first) into two smaller ones across a
reducible point, whenever a valid split point exists; used to normalize a chain of
dominant permutations into minimal reducible pieces.

<a id="schubmult.combinatorics.permutation.s"></a>

#### s

```python
@cache
def s(i)
```

The simple reflection swapping 1-indexed positions ``i`` and ``i + 1``.

<a id="schubmult.combinatorics.pipe_dream"></a>

# schubmult.combinatorics.pipe\_dream

Classical pipe dreams (rectangular grid of crosses/bumps encoding a reduced word),
with conversion to/from `RCGraph`/`WCGraph` and grid symmetries (transpose-like duals,
inversion, vertical reflection).

<a id="schubmult.combinatorics.pipe_dream.PipeDream"></a>

## PipeDream Objects

```python
class PipeDream(PlanarHistory, GridPrint)
```

A pipe dream: pipes enter from the left and exit at the top of a triangular grid of
``CROSS``/``BUMP`` tiles. ``perm`` recovers the induced permutation; ``perm_word`` its
associated (not necessarily reduced) word.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.perm_word"></a>

#### perm\_word

```python
@property
def perm_word()
```

The word read off the crosses in native (left-to-right, top-to-bottom) orientation.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.is_reduced"></a>

#### is\_reduced

```python
@property
def is_reduced()
```

Whether ``perm_word`` is a reduced word for ``perm`` (length equals the number of crosses).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.perm"></a>

#### perm

```python
@property
def perm()
```

The permutation induced by this pipe dream (0-Hecke product of the cross word).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph()
```

Convert to an `RCGraph` (crosses of row ``i`` become that row's column labels).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.to_wc_graph"></a>

#### to\_wc\_graph

```python
def to_wc_graph()
```

Convert to a `WCGraph`, analogous to ``to_rc_graph``.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Build the pipe dream whose crosses are exactly the RC graph's marked positions.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.from_wc_graph"></a>

#### from\_wc\_graph

```python
@classmethod
def from_wc_graph(cls, wc_graph)
```

Build the pipe dream whose crosses are exactly the WC graph's marked positions.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.co_pipe_dream"></a>

#### co\_pipe\_dream

```python
def co_pipe_dream()
```

Swap crosses and bumps under the anti-diagonal reflection (the "co" dual pipe dream).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.co_stinkbat_pipe_dream"></a>

#### co\_stinkbat\_pipe\_dream

```python
def co_stinkbat_pipe_dream()
```

Like ``co_pipe_dream`` but preserving (rather than swapping) cross/bump identity under the reflection.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.inverse_pipe_dream"></a>

#### inverse\_pipe\_dream

```python
def inverse_pipe_dream()
```

Pipe dream for ``~self.perm``, obtained by reflecting bumps/crosses through the anti-diagonal.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.reflect_vertically"></a>

#### reflect\_vertically

```python
def reflect_vertically()
```

Flip the grid top-to-bottom.

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

<a id="schubmult.combinatorics.plactic.Plactic.evacuation"></a>

#### evacuation

```python
def evacuation(n)
```

Schutzenberger evacuation within the alphabet ``1..n``: insert the complemented row word ``n + 1 - a``.

<a id="schubmult.combinatorics.plactic.Plactic.iter_boxes"></a>

#### iter\_boxes

```python
@property
def iter_boxes()
```

Filled cells in row-word order (bottom row to top, left to right).

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

<a id="schubmult.combinatorics.plactic.Plactic.iter_outer_corners"></a>

#### iter\_outer\_corners

```python
@property
def iter_outer_corners()
```

Empty cells that can receive a reverse JDT slide.

<a id="schubmult.combinatorics.plactic.Plactic.iter_inner_corners"></a>

#### iter\_inner\_corners

```python
@property
def iter_inner_corners()
```

Holes of the inner shape that can receive a forward JDT slide.

<a id="schubmult.combinatorics.plactic.Plactic.__init__"></a>

#### \_\_init\_\_

```python
def __init__(word=(), inner_shape=None)
```

Build from rows (tuple of tuples), a prebuilt grid, or empty; ``inner_shape`` gives the skew
holes. The grid is stored with one extra border row and column so outer corners always exist.

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

<a id="schubmult.combinatorics.plactic.Plactic.column_word"></a>

#### column\_word

```python
@property
def column_word()
```

Entries read column by column, each column bottom to top.

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

<a id="schubmult.combinatorics.plactic.Plactic.shape"></a>

#### shape

```python
@property
def shape()
```

Row lengths (filled cells per nonempty row).

<a id="schubmult.combinatorics.plactic.Plactic.skew_shape"></a>

#### skew\_shape

```python
@property
def skew_shape()
```

Return the skew shape as a tuple of (row_length, left_offset) pairs.

<a id="schubmult.combinatorics.plactic.Plactic.from_word"></a>

#### from\_word

```python
@classmethod
def from_word(cls, word)
```

RS insertion tableau of a word.

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

<a id="schubmult.combinatorics.plactic.Plactic.rectify"></a>

#### rectify

```python
def rectify()
```

Jeu de taquin rectification of a skew tableau to a straight shape.

<a id="schubmult.combinatorics.plactic.Plactic.superstandard"></a>

#### superstandard

```python
@classmethod
def superstandard(cls, shape)
```

The standard tableau of the given shape filled ``1, 2, ...`` row by row, left to right.

<a id="schubmult.combinatorics.plactic.Plactic.is_semistandard"></a>

#### is\_semistandard

```python
@property
def is_semistandard()
```

Rows weakly increasing and columns strictly increasing (skipping holes).

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

<a id="schubmult.combinatorics.rc_graph"></a>

# schubmult.combinatorics.rc\_graph

`RCGraph`: reduced (pipe-dream-like) compatible graphs, the central combinatorial model for
Schubert polynomials in this package.

An RC graph is a tuple of rows, row ``i`` (0-indexed) a strictly decreasing tuple of column
labels ``>= i + 1``; reading the rows top-to-bottom, right-to-left within a row, gives a reduced
word for the graph's permutation (``perm``). RC graphs support the crystal structure
(`CrystalGraph`), the KOH/Monk product (`product`), conversions to `PipeDream`/`WCGraph`, and
the squash decomposition used to peel a Grassmannian factor off a general RC graph
(``squash_decomp``/``left_squash``, via `AntiRCGraph`).

<a id="schubmult.combinatorics.rc_graph.RCGraph"></a>

## RCGraph Objects

```python
class RCGraph(WCGraph, CrystalGraph)
```

A reduced compatible graph: a tuple of rows (row ``i`` a strictly decreasing tuple of
column labels ``>= i + 1``) whose concatenated reading order gives a reduced word for `perm`.
Construct via the ``WCGraph``/tuple-of-rows constructor, ``RCGraph.principal_rc(perm, n)``, or
``RCGraph.all_rc_graphs(perm, n)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_elem_sym"></a>

#### is\_elem\_sym

```python
@property
def is_elem_sym()
```

Whether ``perm`` is the identity or a single elementary-symmetric-type permutation
(one descent, code entries all 0 or 1).

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_squash"></a>

#### left\_squash

```python
def left_squash(other_rc)
```

Squash product with ``other_rc`` stacked to the left (via `AntiRCGraph.squash_product`).

<a id="schubmult.combinatorics.rc_graph.RCGraph.squash_decomp"></a>

#### squash\_decomp

```python
@cache
def squash_decomp()
```

Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass.

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_squash_decomp"></a>

#### left\_squash\_decomp

```python
def left_squash_decomp()
```

Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass (left variant of ``squash_decomp``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_full_grass"></a>

#### is\_full\_grass

```python
@property
def is_full_grass()
```

Whether ``perm`` is the identity or Grassmannian with descent at the last row.

<a id="schubmult.combinatorics.rc_graph.RCGraph.args"></a>

#### args

```python
@property
def args() -> tuple
```

Return args for sympy compatibility - prevents traversal into tuple contents.

<a id="schubmult.combinatorics.rc_graph.RCGraph.loc_of_inversion"></a>

#### loc\_of\_inversion

```python
def loc_of_inversion(a, b)
```

``(row, col)`` where the left-to-right inversion ``(a, b)`` sits, via a lookup over all inversions.

<a id="schubmult.combinatorics.rc_graph.RCGraph.index_of_inversion"></a>

#### index\_of\_inversion

```python
def index_of_inversion(a, b)
```

Position in `perm_word` of the left-to-right inversion ``(a, b)``, or ``-1`` if absent.

<a id="schubmult.combinatorics.rc_graph.RCGraph.as_reduced_compatible"></a>

#### as\_reduced\_compatible

```python
def as_reduced_compatible()
```

Return ``(perm_word, compatible_sequence)`` where the sequence entry is each letter's row.

<a id="schubmult.combinatorics.rc_graph.RCGraph.little_bump_desc"></a>

#### little\_bump\_desc

```python
def little_bump_desc()
```

Bump the RC graph at its last descent by one (normalizing first if the graph is too short).

<a id="schubmult.combinatorics.rc_graph.RCGraph.inversions"></a>

#### inversions

```python
def inversions()
```

All left-to-right inversion roots, in reading order.

<a id="schubmult.combinatorics.rc_graph.RCGraph.vex"></a>

#### vex

```python
@property
def vex()
```

A vexillary representative in the same crystal component, obtained by iteratively
extending/zeroing rows until ``perm`` avoids the pattern ``2143``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.hw_tab_rep"></a>

#### hw\_tab\_rep

```python
def hw_tab_rep()
```

``(highest_weight_rc, tableau)`` where ``tableau`` is the Yamanouchi tableau of the highest
weight's shape, reverse-raised back to ``self``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.hw_grass_rep"></a>

#### hw\_grass\_rep

```python
def hw_grass_rep()
```

``(highest_weight_rc, grass)`` pairing the highest weight element with `grass`.

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_chute_moves"></a>

#### all\_chute\_moves

```python
def all_chute_moves()
```

All valid chute moves ``(start, end)`` available on this RC graph (see `ChuteMoveElement`).

<a id="schubmult.combinatorics.rc_graph.RCGraph.chute_lower"></a>

#### chute\_lower

```python
def chute_lower(row_num)
```

Apply a chute move lowering the marked element from ``row_num`` into ``row_num - 1``, or ``None`` if invalid.

<a id="schubmult.combinatorics.rc_graph.RCGraph.chute_raise"></a>

#### chute\_raise

```python
def chute_raise(row_num)
```

Apply a chute move raising the marked element from ``row_num`` into ``row_num + 1``, or ``None`` if invalid.

<a id="schubmult.combinatorics.rc_graph.RCGraph.to_top_rc"></a>

#### to\_top\_rc

```python
def to_top_rc()
```

Push every element as far up (chute-raise) as possible; returns ``(rc, raise_seq)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.to_bottom_rc"></a>

#### to\_bottom\_rc

```python
def to_bottom_rc()
```

Push every element as far down (chute-lower) as possible; returns ``(rc, raise_seq)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_inverse_chute_moves"></a>

#### all\_inverse\_chute\_moves

```python
def all_inverse_chute_moves()
```

All valid inverse chute moves ``(start, end)`` (the reverse direction of ``all_chute_moves``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.to_lowest_weight_demaz"></a>

#### to\_lowest\_weight\_demaz

```python
def to_lowest_weight_demaz()
```

Lowest-weight element reached from ``self``'s highest weight by repeated Demazure lowering; ``(rc, raise_seq)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.extract_demazure_atom"></a>

#### extract\_demazure\_atom

```python
def extract_demazure_atom()
```

Extract the Demazure atom associated with this RC graph.

The Demazure atom is indexed by the right key of the weight tableau,
computed using the Willis algorithm with Earliest Weakly Increasing Subsequences (EWIS).

It consists of all RC graphs with the same permutation and length whose
weight tableaux have the same right key.

Returns a list of RC graphs forming the Demazure atom.

<a id="schubmult.combinatorics.rc_graph.RCGraph.raise_seq_word"></a>

#### raise\_seq\_word

```python
@staticmethod
def raise_seq_word(raise_seq)
```

Collapse consecutive repeated entries out of a raising sequence.

<a id="schubmult.combinatorics.rc_graph.RCGraph.grass"></a>

#### grass

```python
@property
def grass()
```

The Grassmannian RC graph with the same weight as ``self``'s highest weight, reverse-raised back to ``self``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.transition"></a>

#### transition

```python
def transition()
```

One step of the Lascoux-Schutzenberger transition: exchange at the last descent and zero out
rows past ``perm``'s trimcode length.

<a id="schubmult.combinatorics.rc_graph.RCGraph.little_bump"></a>

#### little\_bump

```python
def little_bump(i=None, j=None)
```

Bump the reduced word at the inversion ``(i, j)`` (default: the last-descent inversion) up by
one letter, repeatedly re-reducing until valid; the permutation is unchanged.

<a id="schubmult.combinatorics.rc_graph.RCGraph.little_bump_down"></a>

#### little\_bump\_down

```python
def little_bump_down(i, j)
```

Dual of ``little_bump``: bump the letter at inversion ``(i, j)`` down instead of up.

<a id="schubmult.combinatorics.rc_graph.RCGraph.from_reduced_compatible"></a>

#### from\_reduced\_compatible

```python
@classmethod
def from_reduced_compatible(cls, word, seq, length=None)
```

Build an RC graph from a reduced word and its compatible sequence (row assignment per letter).

<a id="schubmult.combinatorics.rc_graph.RCGraph.from_wc_graph"></a>

#### from\_wc\_graph

```python
@classmethod
def from_wc_graph(cls, wc_graph: WCGraph) -> RCGraph
```

Reduce a `WCGraph` down to an `RCGraph` by canceling matching positive/negative inversion pairs.

<a id="schubmult.combinatorics.rc_graph.RCGraph.crystal_weight"></a>

#### crystal\_weight

```python
@cached_property
def crystal_weight()
```

Alias for ``length_vector``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.tableau_decomp"></a>

#### tableau\_decomp

```python
def tableau_decomp() -> tuple[NilPlactic, Plactic]
```

Split ``self`` into a tuple of RC graphs, one per column of descents (vertical cuts at each descent).

<a id="schubmult.combinatorics.rc_graph.RCGraph.sorted_length_vector"></a>

#### sorted\_length\_vector

```python
@cached_property
def sorted_length_vector()
```

``length_vector`` sorted into weakly decreasing order.

<a id="schubmult.combinatorics.rc_graph.RCGraph.extremal_weight"></a>

#### extremal\_weight

```python
@property
def extremal_weight()
```

The extremal weight of ``self``'s crystal component (via ``_extremal_weight``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.forest_invariant"></a>

#### forest\_invariant

```python
@property
def forest_invariant()
```

The forest (indexed by ``omega_invariant``) attached to ``self`` under omega-insertion.

<a id="schubmult.combinatorics.rc_graph.RCGraph.omega_invariant"></a>

#### omega\_invariant

```python
@property
@cache
def omega_invariant()
```

Omega-insertion of the reversed ``perm_word`` (P-symbol data used for forest/K-theoretic invariants).

<a id="schubmult.combinatorics.rc_graph.RCGraph.w0_automorphism"></a>

#### w0\_automorphism

```python
def w0_automorphism(n=None)
```

Conjugate the reduced word by ``w0`` on ``n`` letters (or the minimal ``n`` fitting ``self``),
producing the RC graph for ``w0 * perm * w0`` on the same crystal-compatible footing.

<a id="schubmult.combinatorics.rc_graph.RCGraph.antiaut"></a>

#### antiaut

```python
def antiaut()
```

Convert to `AntiRCGraph` (the anti-orientation view).

<a id="schubmult.combinatorics.rc_graph.RCGraph.forest_weight"></a>

#### forest\_weight

```python
@cached_property
def forest_weight()
```

``forest_invariant``'s composition, padded to ``len(self)`` (see `schubmult.combinatorics.indexed_forests`).

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_extremal"></a>

#### is\_extremal

```python
@property
def is_extremal() -> bool
```

Whether ``self`` is the (unique) extremal element of its Demazure crystal weight class:
weakly decreasing length vector matching the highest weight, minimal among ties by sorting-permutation length.

<a id="schubmult.combinatorics.rc_graph.RCGraph.demazure_weight"></a>

#### demazure\_weight

```python
@property
def demazure_weight() -> tuple[int, ...]
```

Weight of the distinguished Demazure extremal element.

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_rc"></a>

#### is\_rc

```python
@property
def is_rc() -> bool
```

Whether every entry of row ``i`` (0-indexed) is ``>= i + 1`` (the basic RC graph shape constraint).

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_valid"></a>

#### is\_valid

```python
@property
def is_valid() -> bool
```

Whether ``perm_word`` is reduced for ``perm`` and every entry respects the row-shape constraint.

<a id="schubmult.combinatorics.rc_graph.RCGraph.shiftup"></a>

#### shiftup

```python
def shiftup(shift: int = 1, check_valid=True) -> RCGraph
```

Add ``shift`` to every entry of every row.

<a id="schubmult.combinatorics.rc_graph.RCGraph.right_root_at"></a>

#### right\_root\_at

```python
@cache
def right_root_at(i: int, j: int) -> tuple[int, int]
```

The positive root at grid position ``(i, j)``, transported to the right by the remaining word.

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_root_at"></a>

#### left\_root\_at

```python
@cache
def left_root_at(i: int, j: int) -> tuple[int, int] | None
```

The positive root at grid position ``(i, j)``, transported by everything to its left/above.

<a id="schubmult.combinatorics.rc_graph.RCGraph.inversion_label"></a>

#### inversion\_label

```python
@cache
def inversion_label(i: int, j: int) -> int
```

Row where the inversion ``(i+1, j+1)`` is crossed (its "label" in the weak/Lehmer order sense).

<a id="schubmult.combinatorics.rc_graph.RCGraph.lehmer_label"></a>

#### lehmer\_label

```python
@cache
def lehmer_label(i: int, j: int) -> int
```

Rank of ``inversion_label(i, j)`` among the inversion labels of the roots ``(i', j)``, ``i' <= i``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.perm_word"></a>

#### perm\_word

```python
@cached_property
def perm_word() -> tuple[int, ...]
```

Concatenation of the rows, top to bottom: a reduced word for ``perm``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.reduced_word"></a>

#### reduced\_word

```python
@property
def reduced_word() -> tuple[int, ...]
```

Alias for ``perm_word``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_dom_perm_yamanouchi"></a>

#### is\_dom\_perm\_yamanouchi

```python
def is_dom_perm_yamanouchi(dom_perm: Permutation, perm: Permutation) -> bool
```

Whether ``self`` (assumed to have permutation ``dom_perm``) matches the highest weight of the
Demazure-crystal tensor factor for the ``dom_perm``-part of the product ``S_self.perm * S_dom_perm``
landing on ``perm``, via matching P/weight tableaux against the principal RC graphs.

<a id="schubmult.combinatorics.rc_graph.RCGraph.shape"></a>

#### shape

```python
@property
def shape() -> tuple[int, ...]
```

Row lengths of the Edelman-Greene P-tableau (``p_tableau``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.__invert__"></a>

#### \_\_invert\_\_

```python
def __invert__() -> RCGraph
```

RC graph for ``~perm``, transposed via toggling every marked cell to its mirror position.

<a id="schubmult.combinatorics.rc_graph.RCGraph.normalize"></a>

#### normalize

```python
def normalize() -> RCGraph
```

Resize to ``perm.max_descent`` rows (drop or extend trailing empty rows to the canonical length).

<a id="schubmult.combinatorics.rc_graph.RCGraph.resize"></a>

#### resize

```python
def resize(new_length: int) -> RCGraph
```

Truncate (via ``rowrange``) or extend to exactly ``new_length`` rows.

<a id="schubmult.combinatorics.rc_graph.RCGraph.edelman_greene"></a>

#### edelman\_greene

```python
def edelman_greene() -> tuple[NilPlactic, Plactic]
```

Edelman-Greene correspondence: insert the inversions (in reverse reading order) to build the
``(P, Q)`` pair of a nilCoxeter tableau and a plactic recording tableau.

<a id="schubmult.combinatorics.rc_graph.RCGraph.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other: object) -> object
```

Multiply as elements of the `RCGraphRing` (delegates to that ring's product for `RCGraph` operands).

<a id="schubmult.combinatorics.rc_graph.RCGraph.asdtype"></a>

#### asdtype

```python
def asdtype(cls: type) -> object
```

Convert to the combinatorial-ring element type ``cls`` (via ``cls.dtype().ring.from_rc_graph``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.as_nil_hecke"></a>

#### as\_nil\_hecke

```python
def as_nil_hecke(x: object, y: object | None = None) -> object
```

Represent as a `NilHeckeRing` element: ``polyvalue(x, y) * R(perm)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.has_element"></a>

#### has\_element

```python
@cache
def has_element(i: int, j: int) -> bool
```

Whether row ``i`` (1-indexed) contains the reflection at column ``j`` (i.e. label ``i + j - 1``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.length_vector"></a>

#### length\_vector

```python
@cached_property
def length_vector() -> tuple[int]
```

Row lengths (the crystal weight vector).

<a id="schubmult.combinatorics.rc_graph.RCGraph.lehmer_partial_leq"></a>

#### lehmer\_partial\_leq

```python
@cache
def lehmer_partial_leq(other: RCGraph) -> bool
```

Whether every root's `lehmer_label` in ``self`` is ``<=`` the corresponding label in ``other``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.rowrange"></a>

#### rowrange

```python
def rowrange(start: int, end: int | None = None) -> RCGraph
```

Rows ``[start, end)`` as a fresh RC graph, entries shifted down by ``start``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x: Sequence[Expr],
              y: Sequence[Expr] | None = None,
              *,
              beta: Expr = None,
              prop_beta: bool = False,
              crystal: bool = False) -> Expr
```

Monomial (``y=None``), double (``y`` given), or beta-deformed Grothendieck contribution of this RC graph.

With ``crystal=True``, sums ``polyvalue`` over the whole crystal component instead of just ``self``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.random_rc_graph"></a>

#### random\_rc\_graph

```python
@classmethod
def random_rc_graph(cls, perm: Permutation, length: int = -1) -> RCGraph
```

A uniformly random RC graph for ``perm`` with the given number of rows.

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_rcs_with_word"></a>

#### all\_rcs\_with\_word

```python
@classmethod
def all_rcs_with_word(cls, perm: Permutation,
                      word: tuple[int, ...]) -> set[RCGraph]
```

All RC graphs for ``perm`` whose ``perm_word`` equals ``word`` exactly.

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_rc_graphs"></a>

#### all\_rc\_graphs

```python
@classmethod
def all_rc_graphs(cls,
                  perm: Permutation,
                  length: int = -1,
                  weight: tuple[int, ...] | None = None,
                  *,
                  check_length=False) -> set[RCGraph]
```

All RC graphs for ``perm`` with ``length`` rows (default: ``len(perm.trimcode)``), optionally
restricted to a given ``weight`` (row-length vector). Recursively built via ``pull_out_var``
on the top variable; results are cached by ``(perm, length)`` / ``(perm, weight)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.extend"></a>

#### extend

```python
def extend(extra_rows: int) -> RCGraph
```

Append ``extra_rows`` empty rows at the bottom.

<a id="schubmult.combinatorics.rc_graph.RCGraph.prepend"></a>

#### prepend

```python
def prepend(extra_rows: int) -> RCGraph
```

Insert ``extra_rows`` empty rows at the top (shifting existing entries up accordingly).

<a id="schubmult.combinatorics.rc_graph.RCGraph.pieri_insert"></a>

#### pieri\_insert

```python
def pieri_insert(descent,
                 rows,
                 return_reflections=False,
                 backwards=True,
                 left=False)
```

Insert one crossing per entry of ``rows`` (grouped by row) at the given ``descent``, rectifying
as needed to stay a valid RC graph; the Pieri-rule building block used by `zero_out_last_row`,
`pull_out_row`, and related transition-formula operations.

<a id="schubmult.combinatorics.rc_graph.RCGraph.weight"></a>

#### weight

```python
@property
def weight() -> tuple[int, ...]
```

Flat weight sequence: row index (1-indexed) repeated once per crossing in that row.

<a id="schubmult.combinatorics.rc_graph.RCGraph.perm"></a>

#### perm

```python
@property
def perm() -> Permutation
```

The permutation induced by this RC graph: the reduced product of its reflections.

<a id="schubmult.combinatorics.rc_graph.RCGraph.hecke_perm"></a>

#### hecke\_perm

```python
@property
def hecke_perm() -> Permutation
```

The 0-Hecke (Demazure) product of the reflections, allowing non-length-additive steps.

<a id="schubmult.combinatorics.rc_graph.RCGraph.multiply_reps"></a>

#### multiply\_reps

```python
@classmethod
def multiply_reps(cls, drep1, drep2)
```

Multiply two dicts of ``{tuple-of-RCGraph: coeff}`` representations by squash-producting their
factors in size order, returning an `RCGraphRing` element.

<a id="schubmult.combinatorics.rc_graph.RCGraph.cem_rep"></a>

#### cem\_rep

```python
@cached_property
def cem_rep()
```

``self``'s coefficient in the complete-elementary-monomial (CEM) basis expansion of its own permutation.

<a id="schubmult.combinatorics.rc_graph.RCGraph.custom_cem_rep"></a>

#### custom\_cem\_rep

```python
@cache
def custom_cem_rep(partition)
```

Like ``cem_rep``, but expanding against an explicit dominant ``partition`` rather than the inferred one.

<a id="schubmult.combinatorics.rc_graph.RCGraph.sem_rep"></a>

#### sem\_rep

```python
@cache
def sem_rep(length=None)
```

``custom_cem_rep`` specialized to the staircase partition ``w0(length).trimcode``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.in_CEM_basis"></a>

#### in\_CEM\_basis

```python
@classmethod
@cache
def in_CEM_basis(
    cls,
    perm: Permutation,
    length: int,
    partition: tuple[int] | None = None
) -> dict[RCGraph, dict[tuple[RCGraph], int]]
```

Expand ``S_perm``'s complete-elementary-monomial (CEM) representation, restricted to the
pieces landing on RC graphs of the given ``length``: ``{rc: {tuple-of-elem-sym-RCGraphs: coeff}}``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.full_CEM"></a>

#### full\_CEM

```python
@classmethod
@cache
def full_CEM(
    cls,
    perm: Permutation,
    length: int,
    partition: tuple[int] | None = None
) -> dict[RCGraph, dict[tuple[RCGraph], int]]
```

Like ``in_CEM_basis`` but expanding against the full staircase-bounded strict dominant permutation.

<a id="schubmult.combinatorics.rc_graph.RCGraph.full_double_elem_sym_squash"></a>

#### full\_double\_elem\_sym\_squash

```python
def full_double_elem_sym_squash(p, yvars, zvars)
```

Sum ``double_elem_sym_squash(p, ...)`` over every elementary-symmetric RC graph of degree ``p``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.snap_qy"></a>

#### snap\_qy

```python
def snap_qy()
```

Merge adjacent rows where possible to reach a quasi-Yamanouchi representative, raising
``ValueError`` if that's not achievable.

<a id="schubmult.combinatorics.rc_graph.RCGraph.double_elem_rep"></a>

#### double\_elem\_rep

```python
def double_elem_rep(yvars, size)
```

Express ``self`` in the double elementary-symmetric basis of `BoundedRCFactorAlgebra`, via its
(assumed unique) `full_CEM` decomposition, recursively correcting for lower-order terms.

<a id="schubmult.combinatorics.rc_graph.RCGraph.double_elem_sym_squash"></a>

#### double\_elem\_sym\_squash

```python
def double_elem_sym_squash(weight, _yvars, _zvars)
```

Squash-product ``self`` with the elementary-symmetric RC graph of the given ``weight``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.elem_sym_rcs"></a>

#### elem\_sym\_rcs

```python
@classmethod
def elem_sym_rcs(cls, p, k, length=None, weight=None) -> set[RCGraph]
```

All RC graphs for the elementary-symmetric permutation ``uncode([0]*(k-p) + [1]*p)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.transpose"></a>

#### transpose

```python
def transpose(length: int | None = None) -> RCGraph
```

RC graph for ``~perm``, built by peeling diagonals off the end of each row.

<a id="schubmult.combinatorics.rc_graph.RCGraph.from_array"></a>

#### from\_array

```python
@classmethod
def from_array(cls, arr, min_length=None) -> RCGraph
```

Build an RC graph from a 2D object array (non-``None``/non-zero cells mark crossings).

<a id="schubmult.combinatorics.rc_graph.RCGraph.one_row"></a>

#### one\_row

```python
@classmethod
def one_row(cls, p: int) -> RCGraph
```

The single-row RC graph ``(p, p-1, ..., 1)`` (for the permutation with one nonzero code entry ``p``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.weak_order_leq"></a>

#### weak\_order\_leq

```python
def weak_order_leq(other: RCGraph) -> bool
```

Whether every root's `lehmer_label`/`inversion_label` in ``self`` is ``<=`` in ``other`` (weak order comparison).

<a id="schubmult.combinatorics.rc_graph.RCGraph.w_key_cache"></a>

#### w\_key\_cache

noqa: RUF012

<a id="schubmult.combinatorics.rc_graph.RCGraph.rc_cache"></a>

#### rc\_cache

noqa: RUF012

<a id="schubmult.combinatorics.rc_graph.RCGraph.toggle_ref_at"></a>

#### toggle\_ref\_at

```python
def toggle_ref_at(i: int, j: int) -> RCGraph
```

Add or remove the crossing at 1-indexed grid position ``(i, j)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.principal_rc_factorization"></a>

#### principal\_rc\_factorization

```python
@classmethod
@cache
def principal_rc_factorization(cls, perm: Permutation) -> tuple[RCGraph]
```

Factor the principal RC graph of ``perm`` into a tuple of elementary-symmetric RC graphs,
one per nonzero code entry, peeled off from the top descent down.

<a id="schubmult.combinatorics.rc_graph.RCGraph.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
@cache
def zero_out_last_row() -> RCGraph
```

Drop the (empty) last row, exchanging descents down via `pieri_insert` so the permutation is preserved.

Core step of the Lascoux-Schutzenberger transition formula.

<a id="schubmult.combinatorics.rc_graph.RCGraph.zero_out_last_column"></a>

#### zero\_out\_last\_column

```python
def zero_out_last_column(width) -> RCGraph
```

Transpose analogue of `zero_out_last_row`: drop the last column down to ``width``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.zero_out_in_place"></a>

#### zero\_out\_in\_place

```python
def zero_out_in_place() -> RCGraph
```

Normalize then repeatedly `zero_out_last_row`, restoring the original row count.

<a id="schubmult.combinatorics.rc_graph.RCGraph.alt_product"></a>

#### alt\_product

```python
def alt_product(other)
```

Alternate product of ``self`` and ``other`` (shifted up), searching the zero-action orbit of
``self`` for representatives that stack validly with ``other``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.crystal_length"></a>

#### crystal\_length

```python
def crystal_length() -> int
```

Number of rows.

<a id="schubmult.combinatorics.rc_graph.RCGraph.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(row: int) -> RCGraph | None
```

Crystal lowering operator ``f_row``: pair letters of row ``row`` with larger unpaired letters
of row ``row + 1``, move the least unpaired letter of ``row`` down if that stays a valid RC graph.

<a id="schubmult.combinatorics.rc_graph.RCGraph.quasi_raising_operator"></a>

#### quasi\_raising\_operator

```python
def quasi_raising_operator(row: int) -> RCGraph | None
```

``raising_operator``, but only if it leaves the reduced word ``perm_word`` unchanged.

<a id="schubmult.combinatorics.rc_graph.RCGraph.quasi_lowering_operator"></a>

#### quasi\_lowering\_operator

```python
def quasi_lowering_operator(row: int) -> RCGraph | None
```

``lowering_operator``, but only if it leaves the reduced word ``perm_word`` unchanged.

<a id="schubmult.combinatorics.rc_graph.RCGraph.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(row: int) -> RCGraph | None
```

Crystal raising operator ``e_row``, dual to ``lowering_operator``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.right_zero_act"></a>

#### right\_zero\_act

```python
def right_zero_act() -> set[RCGraph]
```

All RC graphs one row longer that reduce back to ``self`` under `zero_out_last_row`
(the covering set used to build the crystal upward).

<a id="schubmult.combinatorics.rc_graph.RCGraph.bisect_left_coords_index"></a>

#### bisect\_left\_coords\_index

```python
@cache
def bisect_left_coords_index(row: int,
                             col: int,
                             lo: int = 0,
                             hi: int | None = None) -> int
```

Binary search over ``perm_word`` positions for the insertion point of grid coordinate ``(row, col)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.exchange_property"></a>

#### exchange\_property

```python
def exchange_property(descent: int,
                      return_row: bool = False,
                      left: bool = False) -> RCGraph | tuple[RCGraph, int]
```

Toggle off the crossing realizing the simple root ``(descent, descent + 1)``, per the exchange
property; optionally also return the row it was found in.

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_to_right_inversion"></a>

#### left\_to\_right\_inversion

```python
@cache
def left_to_right_inversion(index: int) -> tuple[int, int]
```

Positive root of the ``index``-th letter of ``perm_word``, transported to the right.

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_to_right_left_inversion"></a>

#### left\_to\_right\_left\_inversion

```python
@cache
def left_to_right_left_inversion(index: int) -> tuple[int, int]
```

Positive root of the ``index``-th letter of ``perm_word``, transported to the left.

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_to_right_inversion_coords"></a>

#### left\_to\_right\_inversion\_coords

```python
@cache
def left_to_right_inversion_coords(index: int) -> tuple[int, int]
```

Grid coordinates ``(row, col)`` of the ``index``-th letter of ``perm_word``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.principal_rc"></a>

#### principal\_rc

```python
@classmethod
def principal_rc(cls, perm: Permutation, length: int | None = None) -> RCGraph
```

The canonical (dominant/staircase-filled) RC graph for ``perm``: row ``i`` is
``(i + code[i], ..., i + 1)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.p_tableau"></a>

#### p\_tableau

```python
@cached_property
def p_tableau() -> NilPlactic
```

Edelman-Greene ``P``-tableau (alias for ``edelman_greene()[0]``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.q_tableau"></a>

#### q\_tableau

```python
@cached_property
def q_tableau() -> Plactic
```

Edelman-Greene ``Q``-tableau, recording tableau (alias for ``edelman_greene()[1]``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.weight_tableau"></a>

#### weight\_tableau

```python
@cached_property
def weight_tableau() -> Plactic
```

Plactic tableau recording the RC graph's weight, via column Edelman-Greene RSK insertion.

<a id="schubmult.combinatorics.rc_graph.RCGraph.monk_insert"></a>

#### monk\_insert

```python
def monk_insert(row)
```

Insert a new crossing at ``row`` via the (equivariant) Monk rule, cascading corrections
upward through earlier rows as needed to stay valid.

<a id="schubmult.combinatorics.rc_graph.RCGraph.huang_bump"></a>

#### huang\_bump

```python
def huang_bump(a, b)
```

Toggle off the inversion ``(a, b)`` and re-insert it one column to the right, rectifying any
resulting invalid crossings (a Huang-style bump used in transition-formula bijections).

<a id="schubmult.combinatorics.rc_graph.RCGraph.product"></a>

#### product

```python
@cache
def product(other: RCGraph) -> dict[RCGraph, int]
```

Compute the product of this RC graph with another.

<a id="schubmult.combinatorics.rc_graph.RCGraph.prod_with_rc"></a>

#### prod\_with\_rc

```python
def prod_with_rc(other: RCGraph) -> dict[RCGraph, int]
```

Deprecated: Use product() instead.

<a id="schubmult.combinatorics.rc_graph.RCGraph.bpd_transpose"></a>

#### bpd\_transpose

```python
def bpd_transpose() -> RCGraph
```

Transpose via the `BPD` model: convert to a bumpless pipe dream, transpose that, and convert back.

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_potential_coproduct"></a>

#### is\_potential\_coproduct

```python
def is_potential_coproduct(rc1: RCGraph, rc2: RCGraph) -> bool
```

Whether ``(rc1, rc2)`` could be the two factors of a coproduct term for ``self``: necessary
conditions on descents, Bruhat order, and row-length additivity, checked recursively on vertical cuts.

<a id="schubmult.combinatorics.rc_graph.RCGraph.ring_act"></a>

#### ring\_act

```python
def ring_act(elem: FreeAlgebraElement) -> dict[RCGraph, Expr]
```

Act on ``self`` by a `FreeAlgebraElement` (converted to the word basis), applying each word's
letters right to left via ``act``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.act"></a>

#### act

```python
def act(p: int) -> set[RCGraph]
```

Act by the single free-algebra word letter ``p``: all RC graphs obtained by prepending a new
top row realizing the length-additive product ``uncode([p]) * perm``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.iterative_act"></a>

#### iterative\_act

```python
def iterative_act(p: int, insert: bool = True) -> set[RCGraph]
```

Iterative implementation of ``act``, building up the new top row one letter at a time.

<a id="schubmult.combinatorics.rc_graph.RCGraph.__ge__"></a>

#### \_\_ge\_\_

```python
def __ge__(other: object) -> bool
```

``not (self < other)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.__gt__"></a>

#### \_\_gt\_\_

```python
def __gt__(other: object) -> bool
```

``not (self <= other)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.inv"></a>

#### inv

```python
@property
def inv() -> int
```

Length of ``perm``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.rows"></a>

#### rows

```python
@property
def rows() -> int
```

Number of rows.

<a id="schubmult.combinatorics.rc_graph.RCGraph.width"></a>

#### width

```python
@property
def width() -> int
```

Alias for ``cols``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.height"></a>

#### height

```python
@property
def height() -> int
```

Alias for ``rows``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.compatible_sequence"></a>

#### compatible\_sequence

```python
@property
def compatible_sequence() -> tuple[int, ...]
```

Row index (1-indexed) repeated once per crossing in that row, in reading order (paired with ``perm_word``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.cols"></a>

#### cols

```python
@property
def cols() -> int
```

Number of columns: ``len(perm) - 1``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.leibniz_rep"></a>

#### leibniz\_rep

```python
def leibniz_rep() -> tuple
```

Represent ``self`` as a tuple of permutations via repeated ``shiftcut``, one per row from the bottom.

<a id="schubmult.combinatorics.rc_graph.RCGraph.classify_demazure_crystal"></a>

#### classify\_demazure\_crystal

```python
def classify_demazure_crystal() -> tuple[tuple[int], Permutation]
```

``(highest_weight, sorting_perm)`` classifying ``self``'s Demazure crystal isomorphism class.

<a id="schubmult.combinatorics.rc_graph.RCGraph.demazure_isomorphism_class"></a>

#### demazure\_isomorphism\_class

```python
@property
def demazure_isomorphism_class() -> tuple[tuple[int], Permutation]
```

Alias for classify_demazure_crystal for explicit API usage.

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_hw_rcs"></a>

#### all\_hw\_rcs

```python
@classmethod
@cache
def all_hw_rcs(cls,
               perm: Permutation,
               length: int,
               weight=None) -> set[RCGraph]
```

All distinct highest-weight elements among the RC graphs for ``perm`` at ``length`` rows.

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_forest_rcs"></a>

#### all\_forest\_rcs

```python
@classmethod
@cache
def all_forest_rcs(cls, comp: tuple[int, ...], weight=None) -> set[RCGraph]
```

All RC graphs whose `forest_weight` equals the composition ``comp`` (over every permutation
appearing in the forest-dual expansion of ``comp``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_key_rcs"></a>

#### all\_key\_rcs

```python
@classmethod
@cache
def all_key_rcs(cls, comp: tuple[int, ...], weight=None) -> set[RCGraph]
```

All RC graphs whose `extremal_weight` equals the composition ``comp`` (over every permutation
appearing in the key-dual expansion of ``comp``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.all_lw_rcs"></a>

#### all\_lw\_rcs

```python
@classmethod
@cache
def all_lw_rcs(cls,
               perm: Permutation,
               length: int,
               weight=None) -> set[RCGraph]
```

All distinct lowest-weight elements among the RC graphs for ``perm`` at ``length`` rows.

<a id="schubmult.combinatorics.rc_graph.RCGraph.shiftcut"></a>

#### shiftcut

```python
def shiftcut() -> RCGraph
```

Drop the bottom row and shift every remaining row down by one, discarding entries that would
become non-positive (companion step of ``leibniz_rep``).

<a id="schubmult.combinatorics.rc_graph.RCGraph.divdiff_desc"></a>

#### divdiff\_desc

```python
def divdiff_desc(desc: int) -> set[RCGraph]
```

All RC graphs reachable from exchanging then lowering at descent ``desc`` (a single divided-difference step).

<a id="schubmult.combinatorics.rc_graph.RCGraph.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(u: Permutation) -> set[RCGraph]
```

Apply the divided-difference operator for each simple reflection in a reduced word of ``u``
(from the top descent down), via repeated ``divdiff_desc``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.last_descent_strip"></a>

#### last\_descent\_strip

```python
def last_descent_strip() -> tuple[int, ...]
```

Peel off crossings at (or past) the last descent via ``exchange_property``; returns ``(rc, strip)``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.is_forest_rc"></a>

#### is\_forest\_rc

```python
@property
def is_forest_rc() -> bool
```

Whether ``forest_weight`` matches ``length_vector`` (i.e. ``self`` already realizes its own forest weight).

<a id="schubmult.combinatorics.rc_graph.RCGraph.pull_out_row"></a>

#### pull\_out\_row

```python
def pull_out_row(row: int, keep_size=False) -> tuple[tuple, RCGraph]
```

Remove all crossings of ``row`` (assumed a descent/empty row), reflowing the rows above via
`pieri_insert`/`_pieri_rectify` to preserve validity.

<a id="schubmult.combinatorics.rc_graph.RCGraph.little_bump_zero"></a>

#### little\_bump\_zero

```python
def little_bump_zero()
```

Bump the empty last-descent row by one, normalizing and recursing until the descent is fully cleared.

<a id="schubmult.combinatorics.rc_graph.RCGraph.dualpieri"></a>

#### dualpieri

```python
def dualpieri(mu: Permutation, w: Permutation) -> set[tuple[tuple, RCGraph]]
```

Dual Pieri expansion (RC graph analogue of `schubmult.mult.positivity.dualpieri`): peels one
column of variables at a time via `divdiff_perm`/`pull_out_row`.

<a id="schubmult.combinatorics.rc_graph.RCGraph.divdiff_act_dict"></a>

#### divdiff\_act\_dict

```python
@staticmethod
def divdiff_act_dict(dct, *s_list) -> dict[RCGraph, Expr]
```

Apply `divdiff_desc` for each simple reflection index in ``s_list`` (right to left) to every
RC graph key of ``dct``, accumulating coefficients.

<a id="schubmult.combinatorics.rc_graph.RCGraph.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key: int | tuple[int, int]) -> tuple[int, ...] | int
```

``self[i]`` -> row ``i``; ``self[i, j]`` -> the crossing label at 0-indexed ``(i, j)`` or ``None``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.__lt__"></a>

#### \_\_lt\_\_

```python
def __lt__(other: object) -> bool
```

Compare first by ``perm.trimcode``, then lexicographically by inversion labels.

<a id="schubmult.combinatorics.rc_graph.RCGraph.__le__"></a>

#### \_\_le\_\_

```python
def __le__(other: object) -> bool
```

``self < other or self == other``.

<a id="schubmult.combinatorics.rc_graph.RCGraph.weight_bump"></a>

#### weight\_bump

```python
def weight_bump() -> RCGraph
```

Extend by one row and shift up by one (a crystal-structure-preserving perturbation, used by
`CrystalGraph.weight_reflection`'s fallback).

<a id="schubmult.combinatorics.rc_graph.RCGraph.inverse_crystal_product"></a>

#### inverse\_crystal\_product

```python
def inverse_crystal_product(other) -> RCGraph
```

Product with ``other`` in `RCGraphRing`, projecting each term to its crystal highest weight.

<a id="schubmult.combinatorics.rc_graph.RCGraph.monk_rc"></a>

#### monk\_rc

```python
@classmethod
def monk_rc(cls, row, descent)
```

The elementary-symmetric RC graph of degree 1 marking ``row`` among ``descent`` rows (used by
the double Monk-rule squash helpers).

<a id="schubmult.combinatorics.rc_graph.RCGraph.forest_poly_value"></a>

#### forest\_poly\_value

```python
def forest_poly_value(x: Sequence[Expr],
                      y: Sequence[Expr] | None = None) -> Expr
```

Double polynomial value using the forest/vine-model column convention (``y`` indexed by inversion column).

<a id="schubmult.combinatorics.root_tableau"></a>

# schubmult.combinatorics.root\_tableau

`RootTableau`: root-labeled tableaux implementing dual Knuth equivalence via JDT
(jeu-de-taquin) slides, with the Edelman-Greene invariant preserved by the crystal operators.

<a id="schubmult.combinatorics.root_tableau.RootTableau"></a>

## RootTableau Objects

```python
class RootTableau(CrystalGraph, GridPrint)
```

A tableau of positive roots (grid cells hold ``(root, recording_letter)`` pairs)
implementing dual Knuth equivalence via up/down JDT slides. ``edelman_greene_invariant``
is preserved by the crystal raising/lowering operators.

<a id="schubmult.combinatorics.root_tableau.RootTableau.edelman_greene_invariant"></a>

#### edelman\_greene\_invariant

```python
@property
def edelman_greene_invariant()
```

The Edelman-Greene insertion tableau's row word of the reduced word (computed on the
``w0``-reversed word and mapped back), as a tuple. Constant on crystal components.

<a id="schubmult.combinatorics.root_tableau.RootTableau.eg_root"></a>

#### eg\_root

```python
def eg_root(index)
```

The ``index``-th right root of the EG-invariant reduced word.

<a id="schubmult.combinatorics.root_tableau.RootTableau.eg_row_word"></a>

#### eg\_row\_word

```python
@property
def eg_row_word()
```

Reduced word recovered from the row word of roots by repeatedly peeling off the last simple root.

<a id="schubmult.combinatorics.root_tableau.RootTableau.shape"></a>

#### shape

```python
@property
def shape()
```

Row lengths of the (possibly skew) tableau, omitting empty rows.

<a id="schubmult.combinatorics.root_tableau.RootTableau.root_insert_rsk"></a>

#### root\_insert\_rsk

```python
@classmethod
def root_insert_rsk(cls, reduced_word, compatible_seq)
```

Build the root tableau of a compatible pair: Edelman-Greene insert the reduced word, then
fill each cell of the recording tableau with ``(right root of that letter, compatible letter)``.

<a id="schubmult.combinatorics.root_tableau.RootTableau.recording_tableau"></a>

#### recording\_tableau

```python
@property
def recording_tableau()
```

Grid of positions (in the reduced word) of each cell's root.

<a id="schubmult.combinatorics.root_tableau.RootTableau.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc: RCGraph)
```

Root tableau of an RC graph (its reduced word with the row-index compatible sequence).

<a id="schubmult.combinatorics.root_tableau.RootTableau.roots_before"></a>

#### roots\_before

```python
def roots_before(row, col)
```

Boxes whose letter precedes that of ``(row, col)`` in the reading order.

<a id="schubmult.combinatorics.root_tableau.RootTableau.perm"></a>

#### perm

```python
@property
def perm()
```

The permutation of the reduced word.

<a id="schubmult.combinatorics.root_tableau.RootTableau.delete_box"></a>

#### delete\_box

```python
def delete_box(box)
```

Remove the letter whose root sits in ``box`` (if that root is a Bruhat descent of ``perm``) and
re-insert the shortened compatible pair; ``None`` if it is not a descent.

<a id="schubmult.combinatorics.root_tableau.RootTableau.rectify"></a>

#### rectify

```python
def rectify(randomized=False)
```

Jeu de taquin rectification: slide into inner corners until none remain (first corner by
default, or a random one).

<a id="schubmult.combinatorics.root_tableau.RootTableau.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(row, col, check=False)
```

Reverse JDT slide into the outer corner ``(row, col)`` (grid grows if needed), moving boxes
down/right and shifting roots accordingly; with ``check`` the EG invariant is verified.

<a id="schubmult.combinatorics.root_tableau.RootTableau.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(row, col, check=False)
```

Perform a downward/rightward jeu-de-taquin slide starting from the given
(row, col) hole (0-indexed). Boxes from below or to the right are moved
into the hole, preferring the smaller recording letter when both exist.
Returns a new RootTableau (does not mutate self).

<a id="schubmult.combinatorics.root_tableau.RootTableau.iter_boxes"></a>

#### iter\_boxes

```python
@property
def iter_boxes()
```

Occupied cells in row-major order.

<a id="schubmult.combinatorics.root_tableau.RootTableau.iter_outer_corners"></a>

#### iter\_outer\_corners

```python
@property
def iter_outer_corners()
```

Empty cells that can receive an `up_jdt_slide`.

<a id="schubmult.combinatorics.root_tableau.RootTableau.iter_inner_corners"></a>

#### iter\_inner\_corners

```python
@property
def iter_inner_corners()
```

Empty cells that can receive a `down_jdt_slide`.

<a id="schubmult.combinatorics.root_tableau.RootTableau.is_valid"></a>

#### is\_valid

```python
@property
def is_valid()
```

Whether the reconstructed RC graph is valid.

<a id="schubmult.combinatorics.root_tableau.RootTableau.reduced_word"></a>

#### reduced\_word

```python
@property
def reduced_word()
```

Reduced word read off the grid (see `_word_from_grid`).

<a id="schubmult.combinatorics.root_tableau.RootTableau.compatible_sequence"></a>

#### compatible\_sequence

```python
@property
def compatible_sequence()
```

Compatible sequence read off the grid alongside the reduced word.

<a id="schubmult.combinatorics.root_tableau.RootTableau.word_grid"></a>

#### word\_grid

```python
@property
def word_grid()
```

Grid of reduced-word letters (one per occupied cell).

<a id="schubmult.combinatorics.root_tableau.RootTableau.grid_word"></a>

#### grid\_word

```python
@property
def grid_word()
```

Letters of `word_grid` in row-word order.

<a id="schubmult.combinatorics.root_tableau.RootTableau.order_grid"></a>

#### order\_grid

```python
@property
def order_grid()
```

Grid giving each cell's position in the reduced word.

<a id="schubmult.combinatorics.root_tableau.RootTableau.letter_at"></a>

#### letter\_at

```python
def letter_at(row, col)
```

Reduced-word letter at ``(row, col)``.

<a id="schubmult.combinatorics.root_tableau.RootTableau.__init__"></a>

#### \_\_init\_\_

```python
def __init__(grid, print_only=False)
```

Wrap an object grid of ``(root, letter)`` cells; unless ``print_only``, verify each cell's root
matches the EG-invariant root of its letter.

<a id="schubmult.combinatorics.root_tableau.RootTableau.weight_tableau"></a>

#### weight\_tableau

```python
@property
def weight_tableau()
```

RS insertion tableau of the row word of compatible letters (the crystal weight tableau).

<a id="schubmult.combinatorics.root_tableau.RootTableau.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

Crystal ``epsilon_index`` of the weight tableau.

<a id="schubmult.combinatorics.root_tableau.RootTableau.eg_grid"></a>

#### eg\_grid

```python
@property
def eg_grid()
```

Grid of EG-invariant root indices.

<a id="schubmult.combinatorics.root_tableau.RootTableau.eg_index_word"></a>

#### eg\_index\_word

```python
@property
def eg_index_word()
```

For each cell in row-word order, the index of its root in the EG-invariant word.

<a id="schubmult.combinatorics.root_tableau.RootTableau.row_word"></a>

#### row\_word

```python
@property
def row_word()
```

Compatible letters read bottom row to top, left to right.

<a id="schubmult.combinatorics.root_tableau.RootTableau.root_row_word"></a>

#### root\_row\_word

```python
@property
def root_row_word()
```

Roots read in row-word order.

<a id="schubmult.combinatorics.root_tableau.RootTableau.rc_graph"></a>

#### rc\_graph

```python
@property
def rc_graph()
```

Reconstruct the RC-graph from the root tableau.

<a id="schubmult.combinatorics.root_tableau.RootTableau.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(i)
```

EG-invariant root of the ``i``-th cell in row-word order.

<a id="schubmult.combinatorics.root_tableau.RootTableau.iter_boxes_row_word_order"></a>

#### iter\_boxes\_row\_word\_order

```python
@property
def iter_boxes_row_word_order()
```

Occupied cells bottom row to top, left to right.

<a id="schubmult.combinatorics.root_tableau.RootTableau.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

Crystal ``e_i``: apply the RC graph raising operator, rebuild the root tableau, and slide it
back into this tableau's shape with `up_jdt_slide`; ``None`` if undefined.

<a id="schubmult.combinatorics.root_tableau.RootTableau.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(row)
```

Crystal ``f_row`` computed directly on RC graph rows ``row`` and ``row+1`` by the bracketing rule,
then rebuilt and slid back into shape; ``None`` if undefined. Asserts the EG invariant.

<a id="schubmult.combinatorics.root_tableau.RootTableau.raising_operator_direct"></a>

#### raising\_operator\_direct

```python
def raising_operator_direct(i)
```

Direct crystal raising operator e_i using EG invariant tracking.

Key insight: EG invariant is preserved, but eg_index_word changes.

<a id="schubmult.combinatorics.root_tableau.RootTableau.lowering_operator_direct"></a>

#### lowering\_operator\_direct

```python
def lowering_operator_direct(i)
```

Direct crystal lowering operator f_i using EG invariant tracking.

<a id="schubmult.combinatorics.schub_poly"></a>

# schubmult.combinatorics.schub\_poly

Backwards-compatible re-export shim; see `schubmult.symbolic.common_polys` for the actual implementations.

<a id="schubmult.combinatorics.schubert_monomial_graph"></a>

# schubmult.combinatorics.schubert\_monomial\_graph

Base class for Schubert monomial graph structures.

Provides a common interface for structures like RCGraph and BPD that represent
monomials in Schubert polynomials via grid-based diagrams.

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph"></a>

## SchubertMonomialGraph Objects

```python
class SchubertMonomialGraph(ABC)
```

Abstract base class for Schubert monomial graph structures.

This class provides a common interface for combinatorial objects that:
- Represent monomials in Schubert polynomials
- Have a grid/matrix-like structure with rows and columns
- Correspond to a permutation
- Have a weight (as a vector or tuple)

Concrete implementations include:
- RCGraph: Reduced-compatible graphs with crystal structure
- BPD: Bumpless pipe dreams

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.perm"></a>

#### perm

```python
@property
@abstractmethod
def perm() -> Permutation
```

Return the permutation associated with this monomial graph.

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.rows"></a>

#### rows

```python
@property
@abstractmethod
def rows() -> int
```

Number of rows in the grid representation.

**Returns**:

  Number of rows

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.cols"></a>

#### cols

```python
@property
@abstractmethod
def cols() -> int
```

Number of columns in the grid representation.

**Returns**:

  Number of columns

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.width"></a>

#### width

```python
@property
def width() -> int
```

Width of the grid (alias for cols).

**Returns**:

  Number of columns

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.height"></a>

#### height

```python
@property
def height() -> int
```

Height of the grid (alias for rows).

**Returns**:

  Number of rows

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.permutation"></a>

#### permutation

```python
@property
def permutation() -> Permutation
```

Alias for perm property.

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.__getitem__"></a>

#### \_\_getitem\_\_

```python
@abstractmethod
def __getitem__(key) -> Any
```

Access elements of the grid.

**Arguments**:

- `key` - Index or tuple of indices
  

**Returns**:

  Element(s) at the specified position

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.normalize"></a>

#### normalize

```python
@abstractmethod
def normalize() -> "SchubertMonomialGraph"
```

Return a normalized version of the monomial graph.

Normalization may involve reordering or simplifying the structure
while preserving its combinatorial properties.

**Returns**:

  Normalized SchubertMonomialGraph object

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.polyvalue"></a>

#### polyvalue

```python
@abstractmethod
def polyvalue(x, y=None, **kwargs) -> Expr
```

Compute the polynomial value represented by this monomial graph.

**Returns**:

  Polynomial representation (e.g., as a SymPy expression)

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.right_zero_act"></a>

#### right\_zero\_act

```python
@abstractmethod
def right_zero_act() -> set["SchubertMonomialGraph"]
```

Compute the right action of a zero (adding a row/column).

**Returns**:

  Set of resulting monomial graphs

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.product"></a>

#### product

```python
@abstractmethod
def product(
        other: "SchubertMonomialGraph") -> dict["SchubertMonomialGraph", int]
```

Compute the product of this monomial graph with another.

This represents the Schubert polynomial multiplication at the monomial level.

**Arguments**:

- `other` - Another monomial graph to multiply with
  

**Returns**:

  Dictionary mapping result monomial graphs to their coefficients

<a id="schubmult.combinatorics.set_valued_tableau"></a>

# schubmult.combinatorics.set\_valued\_tableau

`SetValuedTableau`: semistandard set-valued tableaux (Grothendieck-polynomial combinatorics),
with a crystal structure realized via `SetWord`/`SetLetter`.

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

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.rows"></a>

#### rows

```python
@property
def rows()
```

Number of rows (one past the maximum row index present).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.cols"></a>

#### cols

```python
@property
def cols()
```

Number of columns (one past the maximum column index present).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key)
```

``self[row, col]`` -> the label tuple at that box, or ``None`` if empty.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.shape"></a>

#### shape

```python
@property
def shape()
```

Row lengths (number of boxes per row), trailing zeros dropped.

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

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.__eq__"></a>

#### \_\_eq\_\_

```python
def __eq__(other)
```

Equal iff the underlying cell dicts match.

<a id="schubmult.combinatorics.set_word"></a>

# schubmult.combinatorics.set\_word

Set-valued crystal words: `SetLetter` (a subset-of-``{1..n}`` letter with a
GL_n-type crystal structure) and `SetWord` (a tensor of such letters), with
conversions to/from `WCGraph`.

<a id="schubmult.combinatorics.set_word.SetLetter"></a>

## SetLetter Objects

```python
class SetLetter(CrystalGraph, frozenset)
```

A set of ints with a sqrt(gl_n) crystal structure.

<a id="schubmult.combinatorics.set_word.SetLetter.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

The ambient rank ``n`` (number of crystal indices).

<a id="schubmult.combinatorics.set_word.SetLetter.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Weight vector: multiplicity of each value ``1..n`` in the set.

<a id="schubmult.combinatorics.set_word.SetLetter.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

``e_i``: move an element from ``i+1`` to ``i`` if that increases the weight at ``i``, else ``None``.

<a id="schubmult.combinatorics.set_word.SetLetter.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i)
```

``f_i``: the inverse move to ``raising_operator``, or ``None`` if undefined.

<a id="schubmult.combinatorics.set_word.SetWord"></a>

## SetWord Objects

```python
class SetWord(CrystalGraphTensor)
```

A tuple of SetLetters with a sqrt(gl_n) crystal structure.

<a id="schubmult.combinatorics.set_word.SetWord.to_wc_graph"></a>

#### to\_wc\_graph

```python
def to_wc_graph(rows)
```

Convert to a `WCGraph` with the given number of rows: column ``j`` gets a reflection
at each row in ``self.factors[j-1]``.

<a id="schubmult.combinatorics.set_word.SetWord.from_wc_graph"></a>

#### from\_wc\_graph

```python
@classmethod
def from_wc_graph(cls, wc)
```

Inverse of ``to_wc_graph``: build a `SetWord` from a `WCGraph`, one `SetLetter` per column.

<a id="schubmult.combinatorics.wc_graph"></a>

# schubmult.combinatorics.wc\_graph

`WCGraph`: word-compatible graphs, the common tuple-of-rows base shared by `RCGraph` and
`AntiRCGraph`, where the induced permutation is the Demazure (0-Hecke) product of the row word
rather than requiring it to be reduced.

<a id="schubmult.combinatorics.wc_graph.pivot_transition"></a>

#### pivot\_transition

```python
@cache
def pivot_transition(perm2, target_d=None)
```

Recursively apply `Permutation.pivot_transition` over every nonempty pivot subset until the
maximal corner drops to ``target_d`` (default: ``perm2.max_descent``); returns the set of results.

<a id="schubmult.combinatorics.wc_graph.WCGraph"></a>

## WCGraph Objects

```python
class WCGraph(SchubertMonomialGraph, CrystalGraph, GridPrint, tuple)
```

Word-compatible graph.

Internal representation matches RCGraph: a tuple of rows where row i (0-indexed)
is a strictly decreasing tuple of positive integers >= i + 1.

For WCGraph, the associated permutation is the Demazure product of the
concatenated row word (1-indexed simple reflections).

<a id="schubmult.combinatorics.wc_graph.WCGraph.args"></a>

#### args

```python
@property
def args() -> tuple
```

Return args for sympy compatibility - prevents traversal into tuple contents.

<a id="schubmult.combinatorics.wc_graph.WCGraph.__eq__"></a>

#### \_\_eq\_\_

```python
def __eq__(other: object) -> bool
```

Equal iff both are `WCGraph` with the same rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.trans_co_pipe"></a>

#### trans\_co\_pipe

```python
def trans_co_pipe()
```

The complementary graph on twice as many rows: mark every position ``(i+j, j)`` that is
empty in ``self``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.perm_word"></a>

#### perm\_word

```python
@cached_property
def perm_word() -> tuple[int, ...]
```

Concatenation of the rows, top to bottom (not necessarily reduced).

<a id="schubmult.combinatorics.wc_graph.WCGraph.perm"></a>

#### perm

```python
@cached_property
def perm() -> Permutation
```

The permutation induced by this graph: the Demazure (0-Hecke) product of `perm_word`.

<a id="schubmult.combinatorics.wc_graph.WCGraph.hecke_perm"></a>

#### hecke\_perm

```python
@property
def hecke_perm() -> Permutation
```

Alias for `perm` (already a Hecke/Demazure product).

<a id="schubmult.combinatorics.wc_graph.WCGraph.is_rc"></a>

#### is\_rc

```python
@property
def is_rc() -> bool
```

Whether every entry of row ``i`` (0-indexed) is ``>= i + 1`` (the basic row-shape constraint).

<a id="schubmult.combinatorics.wc_graph.WCGraph.is_reduced"></a>

#### is\_reduced

```python
@property
def is_reduced()
```

Whether ``perm_word`` has exactly ``perm.inv`` letters (no Hecke-cancelling excess).

<a id="schubmult.combinatorics.wc_graph.WCGraph.is_valid"></a>

#### is\_valid

```python
@property
def is_valid() -> bool
```

Whether every row is strictly decreasing, respects the row-shape constraint, and the
compatible sequence/word pair is compatible.

<a id="schubmult.combinatorics.wc_graph.WCGraph.shiftup"></a>

#### shiftup

```python
def shiftup(shift: int = 1, check_valid=True) -> WCGraph
```

Add ``shift`` to every entry of every row.

<a id="schubmult.combinatorics.wc_graph.WCGraph.normalize"></a>

#### normalize

```python
def normalize() -> WCGraph
```

Resize to ``perm.max_descent`` rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.resize"></a>

#### resize

```python
def resize(new_length: int) -> WCGraph
```

Truncate (via ``rowrange``) or extend to exactly ``new_length`` rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.rowrange"></a>

#### rowrange

```python
def rowrange(start: int, end: int | None = None) -> WCGraph
```

Rows ``[start, end)`` as a fresh graph, entries shifted down by ``start``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.extend"></a>

#### extend

```python
def extend(extra_rows: int) -> WCGraph
```

Append ``extra_rows`` empty rows at the bottom.

<a id="schubmult.combinatorics.wc_graph.WCGraph.toggle_ref_at"></a>

#### toggle\_ref\_at

```python
def toggle_ref_at(i: int, j: int) -> WCGraph
```

Add or remove the reflection at 1-indexed grid position ``(i, j)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.has_element"></a>

#### has\_element

```python
@cache
def has_element(i: int, j: int) -> bool
```

Whether row ``i`` (1-indexed) contains the reflection at column ``j`` (label ``i + j - 1``).

<a id="schubmult.combinatorics.wc_graph.WCGraph.length_vector"></a>

#### length\_vector

```python
@cached_property
def length_vector() -> tuple[int, ...]
```

Row lengths.

<a id="schubmult.combinatorics.wc_graph.WCGraph.weight"></a>

#### weight

```python
@cached_property
def weight() -> tuple[int, ...]
```

Flat weight sequence: row index (1-indexed) repeated once per reflection in that row.

<a id="schubmult.combinatorics.wc_graph.WCGraph.rows"></a>

#### rows

```python
@property
def rows() -> int
```

Number of rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.cols"></a>

#### cols

```python
@property
def cols() -> int
```

Number of columns: ``len(perm) - 1``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.width"></a>

#### width

```python
@property
def width() -> int
```

Alias for ``cols``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.height"></a>

#### height

```python
@property
def height() -> int
```

Alias for ``rows``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.compatible_sequence"></a>

#### compatible\_sequence

```python
@cached_property
def compatible_sequence() -> tuple[int, ...]
```

Row index (1-indexed) repeated once per reflection in that row, in reading order.

<a id="schubmult.combinatorics.wc_graph.WCGraph.to_mbpd"></a>

#### to\_mbpd

```python
@cache
def to_mbpd(n: int | None = None)
```

The marked bumpless pipedream ``Psi(RCP(self))`` (paper
``writing/mbpd.solve.tex``, Theorem "T: main").

This composes the trivial ``WCGraph -> RCP`` repackaging with the
row-unpop bijection ``Psi``.  ``n`` is the ambient grid size (defaults
to ``len(self.perm) - 1`` padded to fit the graph).

Cached: WCGraphs are immutable and hashable, so the (self, n) round
trip is memoized.

<a id="schubmult.combinatorics.wc_graph.WCGraph.from_mbpd"></a>

#### from\_mbpd

```python
@classmethod
@cache
def from_mbpd(cls, mbpd) -> WCGraph
```

Inverse of :meth:`to_mbpd`: the WCGraph ``RCP(Phi(mbpd))`` obtained
from the row-pop bijection ``Phi`` (paper Theorem "T: main").

Cached on the (hashable) ``mbpd`` argument.

<a id="schubmult.combinatorics.wc_graph.WCGraph.is_principal"></a>

#### is\_principal

```python
@property
def is_principal() -> bool
```

Whether ``perm`` equals ``uncode(length_vector)`` (this graph realizes its own dominant weight).

<a id="schubmult.combinatorics.wc_graph.WCGraph.grove_wcs"></a>

#### grove\_wcs

```python
@classmethod
def grove_wcs(cls, comp, length=None, base_rc=None) -> set[WCGraph]
```

Grove polynomial of ``comp`` built by inverting the omega insertion.

We enumerate the compatible set-valued labelings ``kappa`` of the indexed
forest ``F = weak_composition_to_indfor(comp)`` (the grove definition), and
realize each labeling as a WCGraph by writing down its (word, compatible
sequence) pair *explicitly* -- the published inverse of the set-valued omega
insertion -- and feeding it to ``WCGraph.from_word_compatible``.

Concretely:

* The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
  labeling ``Q`` of the principal reduced RC graph form the omega pair for
  ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
  reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
  so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
  word position ``len(W) - Q(v)`` and carries the reduced letter
  ``ell(v) = W[len(W) - Q(v)]``.
* A labeling ``kappa`` places the letter ``ell(v)`` into every row
  ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
  ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
  a word whose rows are strictly decreasing, i.e. exactly the data
  ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
  values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
  equals ``x^kappa``.

Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
the degree ``-1`` homogenizer recording extra labels beyond one per node.

<a id="schubmult.combinatorics.wc_graph.WCGraph.forest_weight"></a>

#### forest\_weight

```python
@cached_property
def forest_weight()
```

``forest_invariant``'s composition, padded to ``len(self)`` (see `schubmult.combinatorics.indexed_forests`).

<a id="schubmult.combinatorics.wc_graph.WCGraph.grove_weight"></a>

#### grove\_weight

```python
@property
def grove_weight()
```

Length vector of ``grove_invariant`` (the base member's forest weight).

<a id="schubmult.combinatorics.wc_graph.WCGraph.grove_invariant"></a>

#### grove\_invariant

```python
@property
def grove_invariant() -> tuple[int, ...]
```

The composition of the grove containing ``self`` -- the inverse of
:meth:`grove_wcs`.

Every WCGraph of a fixed permutation lies in exactly one grove.  The grove
is generated by the set-valued labelings of a single indexed forest ``F``,
and its *base* member is the distinguished graph whose ``forest_weight``
equals its ``length_vector`` (the minimal, one-element-per-node labeling).
The grove weight is that base's ``forest_weight`` -- i.e. ``F.code`` padded
to ``len(self)``.

For a base member (``forest_weight == length_vector``) the grove weight
coincides with the forest weight; for set-valued members it recovers the
composition of the base rather than the (larger, occurrence-inflated)
forest weight of ``self`` itself.

<a id="schubmult.combinatorics.wc_graph.WCGraph.one_row"></a>

#### one\_row

```python
@classmethod
def one_row(cls, a: int) -> WCGraph
```

The single-row graph ``(a, a-1, ..., 1)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.omega_invariant"></a>

#### omega\_invariant

```python
@property
@cache
def omega_invariant()
```

Omega-insertion of the reversed ``perm_word`` (P-symbol data used for forest/K-theoretic invariants).

<a id="schubmult.combinatorics.wc_graph.WCGraph.forest_invariant"></a>

#### forest\_invariant

```python
@property
def forest_invariant()
```

The forest attached to ``self`` under omega-insertion (first component of ``omega_invariant``).

<a id="schubmult.combinatorics.wc_graph.WCGraph.flipped_co_wc"></a>

#### flipped\_co\_wc

```python
def flipped_co_wc()
```

Reflect crossings through the anti-diagonal one row shorter than ``perm``'s length.

<a id="schubmult.combinatorics.wc_graph.WCGraph.to_reduced_compatible_set_sequence"></a>

#### to\_reduced\_compatible\_set\_sequence

```python
def to_reduced_compatible_set_sequence()
```

Reduce ``perm_word`` to a genuine reduced word, grouping the compatible sequence values that
collapse onto each surviving root into label sets; returns ``(word, set_seq)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.from_reduced_compatible_set_sequence"></a>

#### from\_reduced\_compatible\_set\_sequence

```python
@classmethod
def from_reduced_compatible_set_sequence(cls, word, set_seq, length=None)
```

Build a graph from a reduced word and a set-valued compatible sequence (via `_from_root_dict`).

<a id="schubmult.combinatorics.wc_graph.WCGraph.from_word_compatible"></a>

#### from\_word\_compatible

```python
@classmethod
def from_word_compatible(cls, word, seq, length=None)
```

Build a graph from a word and its (weakly increasing) compatible sequence, validating compatibility.

<a id="schubmult.combinatorics.wc_graph.WCGraph.left_to_right_inversion_coords"></a>

#### left\_to\_right\_inversion\_coords

```python
def left_to_right_inversion_coords(index: int) -> tuple[int, int]
```

Grid coordinates ``(row, col)`` of the ``index``-th letter of ``perm_word``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.left_to_right_inversion"></a>

#### left\_to\_right\_inversion

```python
def left_to_right_inversion(index: int) -> tuple[int, int]
```

Positive root of the ``index``-th letter of ``perm_word``, transported to the right.

<a id="schubmult.combinatorics.wc_graph.WCGraph.left_to_right_hecke_inversion"></a>

#### left\_to\_right\_hecke\_inversion

```python
def left_to_right_hecke_inversion(index: int) -> tuple[int, int]
```

Like `left_to_right_inversion`, transported via the 0-Hecke (Demazure) product instead.

<a id="schubmult.combinatorics.wc_graph.WCGraph.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(i: int, j: int) -> tuple[int, int]
```

The positive root at grid position ``(i, j)``, transported to the right by the remaining word.

<a id="schubmult.combinatorics.wc_graph.WCGraph.right_hecke_root_at"></a>

#### right\_hecke\_root\_at

```python
def right_hecke_root_at(i: int, j: int) -> tuple[int, int]
```

Like `right_root_at`, transported via the 0-Hecke (Demazure) product instead.

<a id="schubmult.combinatorics.wc_graph.WCGraph.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x: Sequence[Expr],
              y: Sequence[Expr] | None = None,
              *,
              beta: Expr = None,
              prop_beta: bool = False,
              crystal: bool = False,
              minus_convention=False) -> Expr
```

Monomial (``y=None``), double (``y`` given), or beta-deformed Grothendieck contribution of this graph.

<a id="schubmult.combinatorics.wc_graph.WCGraph.is_elem_sym"></a>

#### is\_elem\_sym

```python
@property
def is_elem_sym() -> bool
```

Whether ``perm``'s trimcode is all zeros then all ones (elementary-symmetric shape).

<a id="schubmult.combinatorics.wc_graph.WCGraph.crystal_length"></a>

#### crystal\_length

```python
def crystal_length() -> int
```

Number of rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.hecke_invariant"></a>

#### hecke\_invariant

```python
@property
def hecke_invariant()
```

Hecke column-insertion RSK pair ``(P, Q)`` for ``(compatible_sequence, reversed(perm_word))``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.strong_hecke_invariant"></a>

#### strong\_hecke\_invariant

```python
@property
def strong_hecke_invariant()
```

K-theoretic rectification of the diagonal-strip increasing tableau.

Lay the ``perm_word`` out along a single anti-diagonal strip, reading
from bottom-left to top-right, as a skew
:class:`~schubmult.combinatorics.increasing_tableau.IncreasingTableau`,
then repeatedly apply the simultaneous inner-corner K-theoretic down
slide until the tableau is rectified (no inner corners remain). The
rectified increasing tableau is the strong Hecke invariant.

<a id="schubmult.combinatorics.wc_graph.WCGraph.elem_sym_wcs"></a>

#### elem\_sym\_wcs

```python
def elem_sym_wcs(p, k, weight=None)
```

All WC graphs for the elementary-symmetric permutation ``uncode([0]*(k-p) + [1]*p)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.to_rc_pieri"></a>

#### to\_rc\_pieri

```python
def to_rc_pieri()
```

Map this WCGraph to an RCGraph via ``_snap_reduced`` + Pieri insertion.

Snap to the underlying reduced RCGraph, then reinsert the missing
(excess) letters row-by-row: the number of letters missing from row
``i + 1`` is ``self.length_vector[i] - reduced.length_vector[i]``, and
those rows (with multiplicity) are Pieri-inserted at ``perm.max_descent``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight() -> tuple[int, ...]
```

Alias for ``length_vector``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.excess"></a>

#### excess

```python
@property
def excess()
```

How far ``perm_word`` overshoots a reduced word: ``len(perm_word) - perm.inv``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i: int) -> WCGraph | None
```

Crystal raising operator ``e_i``, via the recording tableau of `hecke_invariant`.

<a id="schubmult.combinatorics.wc_graph.WCGraph.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i: int) -> WCGraph | None
```

Crystal lowering operator ``f_i``, dual to `raising_operator`.

<a id="schubmult.combinatorics.wc_graph.WCGraph.sorted_length_vector"></a>

#### sorted\_length\_vector

```python
@cached_property
def sorted_length_vector()
```

``length_vector`` sorted into weakly decreasing order.

<a id="schubmult.combinatorics.wc_graph.WCGraph.extremal_weight"></a>

#### extremal\_weight

```python
@property
def extremal_weight()
```

The extremal weight of ``self``'s strong-Hecke-invariant class, padded to ``len(self)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.groth_to_schub"></a>

#### groth\_to\_schub

```python
@classmethod
@cache
def groth_to_schub(cls, groth_perm: Permutation, beta: Expr)
```

Expand the Grothendieck class of ``groth_perm`` in the (beta-deformed) Schubert basis, via
co-pipe-dreams of every WC graph for ``groth_perm``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.schub_to_groth"></a>

#### schub\_to\_groth

```python
@classmethod
@cache
def schub_to_groth(cls, schub_perm: Permutation, beta: Expr)
```

Expand the Schubert class of ``schub_perm`` in the (beta-deformed) Grothendieck basis, via
co-pipe-dreams of every RC graph for ``schub_perm``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.grove_to_forest"></a>

#### grove\_to\_forest

```python
@classmethod
@cache
def grove_to_forest(cls, comp, beta: Expr)
```

Expand the grove class of composition ``comp`` in the forest-weight basis, via co-pipe-dreams
of every WC graph in the grove.

<a id="schubmult.combinatorics.wc_graph.WCGraph.bisect_left_coords_index"></a>

#### bisect\_left\_coords\_index

```python
@cache
def bisect_left_coords_index(row: int,
                             col: int,
                             lo: int = 0,
                             hi: int | None = None) -> int
```

Binary search over ``perm_word`` positions for the insertion point of grid coordinate ``(row, col)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.vertical_cut"></a>

#### vertical\_cut

```python
def vertical_cut(row: int) -> tuple[WCGraph, WCGraph]
```

Split at ``row`` into two graphs: ``(front, back)`` with ``front`` zeroed down to ``row`` rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.disjoint_union"></a>

#### disjoint\_union

```python
def disjoint_union(rc: WCGraph) -> WCGraph
```

Stack ``rc`` beside ``self`` (shifted so their reflections don't collide), same row count.

<a id="schubmult.combinatorics.wc_graph.WCGraph.is_quasi_yamanouchi"></a>

#### is\_quasi\_yamanouchi

```python
@property
def is_quasi_yamanouchi() -> bool
```

Whether no row can be merged into the row above it (a normal-form condition for `dst`).

<a id="schubmult.combinatorics.wc_graph.WCGraph.dst"></a>

#### dst

```python
@property
def dst()
```

Merge mergeable adjacent rows until reaching a quasi-Yamanouchi ("dominant sorting"-normalized) form.

<a id="schubmult.combinatorics.wc_graph.WCGraph.squash_product"></a>

#### squash\_product

```python
@cache
def squash_product(rc: WCGraph) -> WCGraph
```

Disjoint-union ``self`` and ``rc``, then repeatedly `zero_out_last_row` back down to ``len(self)`` rows.

<a id="schubmult.combinatorics.wc_graph.WCGraph.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
@cache
def zero_out_last_row() -> WCGraph
```

Zero out the (empty) last row, realizing the Weigandt/Lascoux
transition (``writing/weigandt_bumpless.tex``, Theorem "transition").

The graph is transported to a marked bumpless pipedream, the underlying
BPD is resized down by one row (dropping the maximal-corner row), and the
result is transported back.  Weight is preserved; the permutation changes
according to the transition.  Unlike the previous Hecke-insertion route,
this is total: it works for non-core-reduced graphs as well.

Cached: the full MBPD round trip is memoized on the (hashable) graph,
so repeated ``squash_product`` calls reuse the result.

<a id="schubmult.combinatorics.wc_graph.WCGraph.right_zero_act"></a>

#### right\_zero\_act

```python
def right_zero_act() -> set[WCGraph]
```

All WC graphs one row longer that reduce back to ``self`` under `zero_out_last_row`
(the covering set used to build the crystal upward).

<a id="schubmult.combinatorics.wc_graph.WCGraph.principal_wc"></a>

#### principal\_wc

```python
@classmethod
def principal_wc(cls, perm, length)
```

The canonical WC graph for ``perm``: the principal RC graph viewed as a `WCGraph`.

<a id="schubmult.combinatorics.wc_graph.WCGraph.product"></a>

#### product

```python
def product(other: SchubertMonomialGraph) -> dict[WCGraph, int]
```

Compute the product of this WC graph with another.

<a id="schubmult.combinatorics.wc_graph.WCGraph.upieri_insert"></a>

#### upieri\_insert

```python
def upieri_insert(descent,
                  rows,
                  return_reflections=False,
                  backwards=True,
                  left=False)
```

Insert one crossing per entry of ``rows`` at the given ``descent``, rectifying as needed;
WCGraph analogue of `RCGraph.pieri_insert`.

<a id="schubmult.combinatorics.wc_graph.WCGraph.pull_out_var_hecke"></a>

#### pull\_out\_var\_hecke

```python
@classmethod
def pull_out_var_hecke(cls, w: Permutation,
                       k: int) -> set[tuple[tuple[int, ...], Permutation]]
```

Hecke/Demazure analogue of pull_out_var(1, w) for fixed first-row size.

Returns all pairs ``(row, wpp)`` with:
- ``row`` a strictly decreasing tuple of simple reflections of size ``k``
- ``wpp`` a permutation such that
  ``Permutation.ref_product(*row) @ wpp.shiftup(1) == w``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.all_wc_graphs_slow"></a>

#### all\_wc\_graphs\_slow

```python
@classmethod
def all_wc_graphs_slow(cls,
                       perm: Permutation,
                       length: int = -1,
                       weight: tuple[int, ...] | None = None) -> set[WCGraph]
```

Generate all WC graphs with Demazure permutation ``perm`` and fixed row count.

The recursion mirrors ``all_rc_graphs`` but uses Hecke/Demazure pull-out on the
first row via ``pull_out_var_hecke``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.all_wc_graphs"></a>

#### all\_wc\_graphs

```python
@classmethod
def all_wc_graphs(cls,
                  perm: Permutation,
                  length: int | None = None,
                  weight: tuple[int, ...] | None = None,
                  *,
                  check_length=False,
                  do_cache=True) -> set[WCGraph]
```

All WC graphs for ``perm`` with ``length`` rows (default: ``len(perm.trimcode)``), optionally
restricted to a given ``weight``. Recursively built via ``pull_out_var`` on the top variable,
enumerating every addable-descent row via a stack-based search; results are cached by
``(perm, length)``/``(perm, weight)`` unless ``do_cache=False``. See also the reference
(slower) implementation `all_wc_graphs_slow`.

<a id="schubmult.combinatorics.wc_graph.WCGraph.grothendieck_polynomial_via_wc"></a>

#### grothendieck\_polynomial\_via\_wc

```python
@classmethod
def grothendieck_polynomial_via_wc(cls,
                                   perm: Permutation,
                                   x: Sequence[Expr],
                                   beta: Expr,
                                   length: int = -1)
```

Compute a Grothendieck candidate by summing WC graph monomials.

Uses weight monomials with ``beta`` exponent ``|word|-inv(perm)``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key: int | tuple[int, int]) -> tuple[int, ...] | int
```

``self[i]`` -> row ``i``; ``self[i, j]`` -> the crossing label at 0-indexed ``(i, j)`` or ``None``.

<a id="schubmult.mult"></a>

# schubmult.mult

Multiplication algorithms for Schubert polynomials.

This module provides kernels for computing products of Schubert polynomials
in various settings:

- schubmult_py: Ordinary (single) Schubert polynomial multiplication
- schubmult_double: Double Schubert polynomial multiplication
- schubmult_q: Quantum Schubert polynomial multiplication
- schubmult_q_double: Quantum double Schubert polynomial multiplication
- grothmult_double: Double Grothendieck multiplication by a degree-one class

Also includes positivity utilities (posify, compute_positive_rep) for root-based representations.

<a id="schubmult.mult._accel"></a>

# schubmult.mult.\_accel

C++ multiplication kernels (the ``schubmult_cpp`` extension, built from ``cpp/`` by setup.py).

``schubmult_py``, ``schubmult_double``, ``schubmult_q_fast``, ``schubmult_q_double_fast`` and the
``*_from_elems`` kernels dispatch here. The extension is a required part of the package; the
pure-Python kernels remain only as the fallback for permutations beyond the compiled MAXN (the
wrappers return ``None`` in that case). Set ``SCHUBMULT_NO_CPP=1`` to force the Python kernels.

<a id="schubmult.mult._accel.schubmult_py"></a>

#### schubmult\_py

```python
def schubmult_py(perm_dict, v)
```

C++ ``schubmult_py``; returns ``None`` if any coefficient is non-integer or the size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_double"></a>

#### schubmult\_double

```python
def schubmult_double(perm_dict, v, var2, var3)
```

C++ ``schubmult_double``; returns ``None`` if the size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_q_fast"></a>

#### schubmult\_q\_fast

```python
def schubmult_q_fast(perm_dict, v, q_var)
```

C++ ``schubmult_q_fast``; non-integer coefficients are handled by linearity, one key at a time.

<a id="schubmult.mult._accel.schubmult_q_double_fast"></a>

#### schubmult\_q\_double\_fast

```python
def schubmult_q_double_fast(perm_dict, v, var2, var3, q_var)
```

C++ ``schubmult_q_double_fast``; returns ``None`` if the size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_double_from_elems"></a>

#### schubmult\_double\_from\_elems

```python
def schubmult_double_from_elems(perm_dict, v, var2, var3, elem_func)
```

C++ ``schubmult_double`` with a custom elementary-symmetric ``elem_func``; ``None`` if size exceeds ``MAXN``.

<a id="schubmult.mult._accel.schubmult_double_alt_from_elems"></a>

#### schubmult\_double\_alt\_from\_elems

```python
def schubmult_double_alt_from_elems(perm_dict, v, var2, var3, elem_func)
```

C++ ``schubmult_double_alt_from_elems``; ``None`` if size exceeds ``MAXN``.

<a id="schubmult.mult.double"></a>

# schubmult.mult.double

Double Schubert polynomial multiplication.

Implements the ``schubmult_double`` kernel: the product of a linear
combination of double Schubert polynomials ``S_u(x, var2)`` with a single
``S_v(x, var3)``, returned as a coefficient dict ``{w: coeff}`` of polynomials
in ``var2``/``var3``. Uses the same ``theta``/``vmu``/v-path recursion as
``schubmult.mult.single`` (see that module), replacing the plain elementary
symmetric contribution with ``elem_sym_func``, which carries the secondary
variables ``var2``/``var3``.

Also provides the "alt"/"from_elems" variants (building the product one
descent-pulled variable at a time via ``pull_out_var``, generic to any choice
of elementary-symmetric-like function), ``nilhecke_mult`` (nilHecke ring
multiplication), and ``schub_coprod_double`` (the double Schubert coproduct).

<a id="schubmult.mult.double.count_sorted"></a>

#### count\_sorted

```python
def count_sorted(mn, tp)
```

Count occurrences of ``tp`` in the sorted sequence ``mn`` via binary search.

<a id="schubmult.mult.double.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var2=None)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the single variable ``x_varnum``.

Equivariant Monk rule: the diagonal term contributes ``var2[u(varnum)]``
(localization of ``x_varnum`` at ``u``) and the off-diagonal terms are the
same Bruhat-cover moves as the ordinary (non-equivariant) ``single_variable``
in ``schubmult.mult.single``.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `varnum` - 1-indexed variable index ``k``.
- `var2` - Secondary (``y``) generating set.
  

**Returns**:

- `dict` - The updated coefficient dict.

<a id="schubmult.mult.double.single_variable_down"></a>

#### single\_variable\_down

```python
def single_variable_down(coeff_dict, varnum, var2=None)
```

Down (descent) variant of ``single_variable``, using ``elem_sym_perms_op``.

<a id="schubmult.mult.double.mult_poly_double"></a>

#### mult\_poly\_double

```python
def mult_poly_double(coeff_dict, poly, var_x=None, var_y=None)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_py`` with
the extra secondary alphabet ``var_y``.

<a id="schubmult.mult.double.mult_poly_double_alt"></a>

#### mult\_poly\_double\_alt

```python
def mult_poly_double_alt(coeff_dict, poly, var_x=None, var_y=None)
```

Variant of ``mult_poly_double`` that folds each factor via ``schubmult_double_dict``
instead of ``single_variable``, so ``poly`` is only ever expanded one variable/factor
at a time in the ``S_v`` basis rather than left as a raw scalar multiplier.

<a id="schubmult.mult.double.mult_poly_down"></a>

#### mult\_poly\_down

```python
def mult_poly_down(coeff_dict, poly)
```

Down (descent) variant of ``mult_poly_double``, using ``single_variable_down``
and the fixed default alphabet ``_vars.var1``.

<a id="schubmult.mult.double.nilhecke_mult"></a>

#### nilhecke\_mult

```python
def nilhecke_mult(coeff_dict1, coeff_dict2)
```

NilHecke ring product of ``coeff_dict1`` (polynomial coefficients) and
``coeff_dict2`` (permutation coefficients acting as divided-difference operators).

For each ``w`` in ``coeff_dict2`` its coefficient polynomial is pushed through
``mult_poly_down`` against ``coeff_dict1``, and each resulting permutation ``v``
is right-multiplied by ``w`` whenever that multiplication is length-additive.

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.double.schubmult_double_pair"></a>

#### schubmult\_double\_pair

```python
@cache
def schubmult_double_pair(perm1, perm2, var2=None, var3=None)
```

``schubmult_double`` specialized to a single ``perm1`` with coefficient 1, cached.

<a id="schubmult.mult.double.schubmult_double_pair_generic"></a>

#### schubmult\_double\_pair\_generic

```python
@cache
def schubmult_double_pair_generic(perm1, perm2)
```

``schubmult_double_pair`` with the fixed generic secondary alphabets ``_vars.var_g1``/``_vars.var_g2``.

<a id="schubmult.mult.double.schubmult_double_pair_generic_alt"></a>

#### schubmult\_double\_pair\_generic\_alt

```python
@cache
def schubmult_double_pair_generic_alt(perm1, perm2)
```

Like ``schubmult_double_pair_generic`` but computed via ``schubmult_double_alt_from_elems``
with the factorial elementary symmetric function, then expanded/simplified.

<a id="schubmult.mult.double.schubmult_double_dict"></a>

#### schubmult\_double\_dict

```python
def schubmult_double_dict(perm_dict1, perm_dict2, var2=None, var3=None)
```

Multiply two coefficient dicts of double Schubert polynomials together.

Computes ``(sum_u coeff1_u S_u(x, var2)) * (sum_v coeff2_v S_v(x, var3))``
by summing ``schubmult_double(perm_dict1, v, var2, var3)`` scaled by
``coeff2_v`` over ``v`` in ``perm_dict2``.

<a id="schubmult.mult.double.schubmult_double"></a>

#### schubmult\_double

```python
def schubmult_double(perm_dict, v, var2=None, var3=None)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the double Schubert polynomial ``S_v(x, var3)``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available (and both
secondary alphabets are given), falling back to the pure-Python
implementation ``_schubmult_double_python`` otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the Schubert polynomial to
  multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}`` (polynomials in ``var2``/``var3``).

<a id="schubmult.mult.double.schubmult_double_alt"></a>

#### schubmult\_double\_alt

```python
def schubmult_double_alt(perm_dict, v, var2=None, var3=None, index=1)
```

Alternate double Schubert product, built by peeling one variable of ``~v`` at a
time via ``pull_out_var`` instead of the ``theta``/v-path recursion.

Multiplies ``sum_u coeff_u S_u(x, var2)`` by ``S_v(x, var3)``, recursing on
``~new_v`` with the elementary symmetric factor coming from
``elem_sym_positional_perms`` at each step.

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_forwards"></a>

#### schubmult\_double\_alt\_from\_elems\_forwards

```python
def schubmult_double_alt_from_elems_forwards(perm_dict,
                                             v,
                                             var2=None,
                                             var3=None,
                                             index=1,
                                             elem_func=None)
```

``schubmult_double_alt`` generalized to an arbitrary elementary-symmetric-like
``elem_func(p, k, x_vars, y_vars)``, processing variables of ``~v`` from the first
pulled-out index forward.

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_backwards"></a>

#### schubmult\_double\_alt\_from\_elems\_backwards

```python
def schubmult_double_alt_from_elems_backwards(perm_dict,
                                              v,
                                              var2=None,
                                              var3=None,
                                              elem_func=None)
```

Like ``schubmult_double_alt_from_elems_forwards`` but processing ``~v``'s pulled-out
variables from the last descent backward, multiplying the elementary-symmetric
factor in *before* recursing (dispatches to the compiled kernel when available).

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_backwards_backwards"></a>

#### schubmult\_double\_alt\_from\_elems\_backwards\_backwards

```python
def schubmult_double_alt_from_elems_backwards_backwards(
        perm_dict, v, var2=None, var3=None, elem_func=None)
```

Variant of ``_schubmult_double_alt_from_elems_backwards_python`` without the
per-``new_v`` memoization cache, recursing on the interim dict instead of the
original ``perm_dict`` at each pulled-out variable.

<a id="schubmult.mult.double.schubmult_double_from_elems"></a>

#### schubmult\_double\_from\_elems

```python
def schubmult_double_from_elems(perm_dict,
                                v,
                                var2=None,
                                var3=None,
                                elem_func=None)
```

``schubmult_double`` generalized to an arbitrary elementary-symmetric-like
``elem_func``, via the ``theta``/v-path recursion (rather than ``pull_out_var``).

Dispatches to the compiled kernel when available, falling back to
``_schubmult_double_from_elems_python``.

<a id="schubmult.mult.double.schubmult_double_down"></a>

#### schubmult\_double\_down

```python
def schubmult_double_down(perm_dict, v, var2=None, var3=None)
```

Down (descent) variant of ``_schubmult_double_python``, using ``elem_sym_perms_op``.

<a id="schubmult.mult.double.schub_coprod_double"></a>

#### schub\_coprod\_double

```python
def schub_coprod_double(mperm, indices, var2=None, var3=None)
```

Coproduct of the double Schubert polynomial ``S_mperm`` restricted to the
variable split named by ``indices``.

Analogue of ``schub_coprod_py``: multiplies the Grassmannian permutation for
``indices`` against ``mperm`` (via ``schubmult_double`` with a merged ``2N``
variable alphabet), splits each resulting permutation's window, and
substitutes the merged alphabet back to ``var2``/``var3``.

**Arguments**:

- `mperm` - Permutation (or array-form list) to take the coproduct of.
- `indices` - Iterable of 1-indexed positions selecting the variable split.
- `var2` - Secondary alphabet for the first factor's variables.
- `var3` - Secondary alphabet for the second factor's variables.
  

**Returns**:

- `dict` - Mapping ``{(firstperm, secondperm): coeff}``.

<a id="schubmult.mult.groth"></a>

# schubmult.mult.groth

beta-Grothendieck Chevalley formula: multiplication of a (single) Grothendieck
polynomial by a bare x_k variable, i.e. x_k * G_w^(beta).

Non-equivariant (y=0) for now; ``GrothendieckRing`` has no coefficient/y genset yet.

Derived from M. Willems, "A Chevalley formula in equivariant K-theory"
(arXiv:math/0603220), Theorem 5 (the ordinary, non-equivariant specialization of
his equivariant Chevalley formula, Theorem 4). Willems indexes K-theory classes
O_w by the *dimension* of the Schubert variety, dual to the *codimension*
indexing used by Schubert/Grothendieck polynomials S_w/G_w; the w0-conjugation
below (``hat_w = w0*w`` going in, ``w0*v`` coming out) translates between the two
conventions. The beta-grading (beta^(d-1) per length difference d = l(v)-l(w))
matches this codebase's beta-deformed Grothendieck polynomial normalization
(beta=0 recovers the classical double Schubert Monk formula). Calibrated against
grothendieck_poly()/to_groth() (see session notes).

<a id="schubmult.mult.groth.chevalley_x_k"></a>

#### chevalley\_x\_k

```python
def chevalley_x_k(w, k, beta, n=None)
```

Coefficients of ``x_k * G_w^(beta)`` in the Grothendieck basis, as a dict
``{v: coeff}`` (``w`` itself never appears: the self-term cancels identically).

<a id="schubmult.mult.groth.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by the single variable ``x_varnum``
(Grothendieck Chevalley formula), via ``chevalley_x_k``.

<a id="schubmult.mult.groth.mult_poly_groth"></a>

#### mult\_poly\_groth

```python
def mult_poly_groth(coeff_dict, poly, var_x, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable_groth``.

<a id="schubmult.mult.groth_double"></a>

# schubmult.mult.groth\_double

K-theoretic Monk formula for double Grothendieck polynomials.

Implements Lenart--Postnikov, *Affine Weyl groups in K-theory and representation
theory* (arXiv:math/0309207), Corollary 8.2 (the :math:`K_T`-Monk formula), in
type :math:`A` and transported to the ``schubmult`` conventions.

Reconciling the conventions
---------------------------
``schubmult`` has no exponentials: it uses the multiplicative formal group law
``a (+) b = a + b + beta*a*b`` with formal inverse ``(-)b = -b/(1 + beta*b)``.
The double Grothendieck polynomial of a simple reflection is

.. math::

    \mathfrak{G}_{s_k}(x, y)
        = \frac{1}{\beta}\Bigl(\prod_{i=1}^{k}(1 + \beta x_i)(1 + \beta y_i) - 1\Bigr),

which at ``beta = -1`` is precisely the class ``1 - x^{w_0(omega_k)} e^{-omega_k}``
of Lemma 8.1(a).  The dictionary between the paper and this package is

* torus characters: ``x^{eps_i}  <->  1 + beta*y_i`` (so that ``(x^{eps_a - eps_b} - 1)/beta``
  is the root ``y_a (-) y_b``, matching ``DoubleGrothendieckRing.exp_root``);
* basis elements: ``[O_{X_{w_0 w}}]  <->  (-beta)^{l(w)} G_w``, since the paper
  indexes structure sheaves by the *dimension* of the Schubert variety whereas
  ``G_w`` has lowest term ``S_w`` of degree ``l(w)`` (codimension indexing).

Under ``u -> w_0 u`` the saturated *decreasing* chains of Corollary 8.2 become
saturated *increasing* chains, the reflections ``t_{ij}`` are unchanged, and the
character prefactor ``x^{nu(J)} = x^{w_0(omega_k) - u(omega_k)}`` (constant in
``J`` because ``omega_k`` is minuscule in type ``A``) becomes
``prod_{i<=k} (1 + beta*y_i)/(1 + beta*y_{u(i)})``.  Rescaling by
``(-beta)^{l(.)}`` turns the signs ``(-1)^{|J|}`` into powers ``beta^{|J|}`` and
yields

.. math::

    \mathfrak{G}_u(x, y)\,\mathfrak{G}_{s_k}(x, z)
        = \frac{1}{\beta}\Bigl(\Theta_u \sum_J \beta^{|J|}\,
          \mathfrak{G}_{u\,r_J}(x, y) - \mathfrak{G}_u(x, y)\Bigr),
    \qquad
    \Theta_u = \prod_{i=1}^{k}\frac{1 + \beta z_i}{1 + \beta y_{u(i)}},

the sum being over the subsets ``J`` of a reduced ``(-omega_k)``-chain of
reflections whose reflections build a saturated increasing Bruhat chain from
``u`` (the empty subset included).  The two secondary alphabets are handled by
``G_{s_k}(x, z) = C G_{s_k}(x, y) + (C - 1)/beta`` with
``C = prod_{i<=k}(1 + beta*z_i)/(1 + beta*y_i)``, which is exactly what turns
the ``y_i`` of ``x^{nu}`` into the ``z_i`` of ``Theta_u``.

By Corollary 15.4 of the same paper a reduced ``(-omega_k)``-chain of
reflections in type ``A_{n-1}`` is

    ``t_{1,n}, t_{1,n-1}, ..., t_{1,k+1}, t_{2,n}, ..., t_{2,k+1}, ..., t_{k,k+1}``.

<a id="schubmult.mult.groth_double.monk_chain"></a>

#### monk\_chain

```python
def monk_chain(k)
```

Reduced ``(-omega_k)``-chain of reflections in ``A_{n-1}``, as ``(i, j)``, ``i <= k < j``.

``omega_k = eps_1 + ... + eps_k``, so this is ``epsilon_chain`` on ``{1, ..., k}``:
every root ``alpha_{ij}`` with ``i, j <= k`` pairs to zero and drops out, and the
survivors all sit at level ``1``.  The order is ``i`` decreasing, then ``j`` decreasing.

<a id="schubmult.mult.groth_double.epsilon_chain"></a>

#### epsilon\_chain

```python
def epsilon_chain(positions, inverse=False, ambient_rank=None)
```

Reduced ``(-eps_A)``-chain of reflections in ``A_{n-1}``, ``A = positions``.

``positions`` is a single index or an iterable of them (repeats allowed), and
``eps_A = sum_{i in A} eps_i``.  With ``inverse=True`` the chain is for ``+eps_A``
instead, which is the weight of the inverse class ``prod_{i in A}(1 + beta*x_i)^{-1}``.

``(-omega_k)``-chains only see the roots ``alpha_{ij}`` with ``i <= k < j``, which is
why ``monk_chain`` multiplies by the whole product ``prod_{i<=k}(1 + beta*x_i)``.
Selecting an arbitrary set of variables needs ``eps_A`` instead, whose chain also
involves the roots ``alpha_{ik}`` with ``i < k``.

Built by Prop. 6.7: the reflections ``s_{alpha, m}`` separating the fundamental
alcove from ``A_{eps_A}``, ordered by the lexicographic key
``(lambda, alpha^vee)^{-1} (-m, (omega_1, alpha^vee), ..., (omega_{n-1}, alpha^vee))``.
Entries are ``(a, b, m)`` for the positive root ``alpha_{ab} = eps_a - eps_b``, ``a < b``;
``m > 0`` means ``b(r) = -alpha`` is negative and the step carries a sign in Thm 6.1.

Concatenating the individual ``(-eps_i)``-chains would also be legal (Prop. 12.2) but
only after translating the blocks, which shifts their levels; going through Prop. 6.7
avoids that and is reduced.  Note a single ``k`` gives ``(i, k, 0)`` for ``i < k`` and
``(k, j, 1)`` for ``j > k``, while for ``A = {1, ..., k}`` every ``alpha_{ij}`` with
``i, j <= k`` pairs to zero and drops out, leaving exactly ``monk_chain(k)``.
Flipping to ``inverse=True`` exchanges those two families.

<a id="schubmult.mult.groth_double.one_plus_beta_x_groth"></a>

#### one\_plus\_beta\_x\_groth

```python
def one_plus_beta_x_groth(coeff_dict,
                          positions,
                          var2=None,
                          beta=None,
                          inverse=False)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by ``prod_{i in positions} (1 + beta*x_i)``.

The one-pass Pieri rule of Theorem 6.1 at ``lambda = -eps_A``; see
``_one_plus_beta_x_terms`` for the coefficient.  ``inverse=True`` gives the inverse
operator ``prod_{i in positions} (1 + beta*x_i)^{-1}``, i.e. ``lambda = +eps_A``.

<a id="schubmult.mult.groth_double.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, var2=None, beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by the single variable ``x_varnum``.

Returns ``{w: coeff_w}``.  This is ``_one_plus_beta_x_terms`` with the
identity subtracted off and ``beta`` divided out; the diagonal coefficient
collapses to ``(-) var2[u(varnum)] = -y/(1 + beta*y)``, the formal inverse of
``var2[u(varnum)]``, which is the localization of ``x_varnum`` at ``u``.

<a id="schubmult.mult.groth_double.elem_sym_perms_groth"></a>

#### elem\_sym\_perms\_groth

```python
def elem_sym_perms_groth(u, k)
```

K-theoretic analogue of ``elem_sym_perms``: ``{w: {d: multiplicity}}``.

Same recursion as ``elem_sym_perms(u, p, k)`` -- a step is any Bruhat cover
``w -> w t_{ij}`` with ``i <= k < j``, and ``j`` is required to weakly decrease along
the chain -- but with the ``p`` cut-off dropped, so chains of every length are
produced and the degree cut-off is left to the coefficient.

This is deliberately *not* a subset-of-a-fixed-chain enumeration.  A
``lambda``-chain imposes a total order on the transpositions, which loses covers such
as ``id < [1,3,2] < [2,3,1]`` (that needs ``t_{23}`` before ``t_{13}``); covers come
from arbitrary upward transpositions, and only ``j`` is constrained.

``d = l(w) - l(u)`` is the chain length.  A position ``i <= k`` may be stepped on more
than once, and distinct chains can land on the same ``w`` at the same ``d``, which is
the source of the K-theoretic multiplicities.

<a id="schubmult.mult.groth_double.grothmult_double_block"></a>

#### grothmult\_double\_block

```python
def grothmult_double_block(coeff_dict,
                           positions,
                           zvar=None,
                           var2=None,
                           beta=None,
                           fgl=True)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by a linear block over ``positions``:

    fgl=True  ->  prod_{i in A} (x_i (+) zvar),   x (+) z = x*(1 + beta*z) + z
    fgl=False ->  prod_{i in A} (x_i - zvar)

``positions`` is an arbitrary index set (repeats allowed), matching the ``index_list``
that ``pull_out_var`` produces, so this is the ``G``-basis analogue of the top-degree
mixed-variable block driving ``schubmult_double_alt`` / ``DoubleSchubertRing.elem_mul``.

``fgl=False`` is the plain ``beta = 0`` block ``(x_1 - z)(x_2 - z)...`` -- still a
perfectly good operator on the ``G`` basis, and the two are interchangeable via
``x (+) z = (1 + beta*z) * (x - (-)z)``, so either can be recovered from the other by
rescaling ``zvar``.

Computed by folding ``single_variable_groth`` one position at a time, using
``(a x_i + b) F = a (x_i F) + b F``.  That keeps every intermediate coefficient
polynomial in ``beta``; the one-pass alternative would expand
``prod_i ((1 + beta*x_i)(1 + beta*z) - 1) / beta**|A|`` by inclusion-exclusion over the
subsets of ``A`` (each term a ``one_plus_beta_x_groth`` call) and only cancel the
``beta^{-|A|}`` at the very end.

<a id="schubmult.mult.groth_double.groth_elem_sym_poly"></a>

#### groth\_elem\_sym\_poly

```python
def groth_elem_sym_poly(p, k, zvar, var_x, beta)
```

``E_p^beta(x_1..x_k; z) = e_p(x_1 (+) z, ..., x_k (+) z)``, ``x (+) z = x(1 + beta*z) + z``.

The double Grothendieck elementary symmetric: ``p == k`` gives
``(x_1(1 + beta*z) + z) ... (x_k(1 + beta*z) + z)`` and ``beta == 0`` gives the
factorial elementary symmetric ``elem_sym_poly(p, k, x, [-z])``.

<a id="schubmult.mult.groth_double.grothmult_double_pieri"></a>

#### grothmult\_double\_pieri

```python
def grothmult_double_pieri(coeff_dict,
                           p,
                           k,
                           zvar=None,
                           var_x=None,
                           var2=None,
                           beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by ``groth_elem_sym_poly(p, k, zvar, var_x, beta)``.

Exact, by folding the ``e_p`` DP over coefficient dicts one shifted variable
``x_i (+) zvar`` at a time (``single_variable_groth`` per step) -- no symbolic
expansion.  A closed-form Pieri rule in the style of ``dom_groth`` -- paths from
``elem_sym_perms_groth`` plus an ``elem_sym_poly`` in the localizations -- is *not*
implemented: grading the paths by ``beta^{d - m}`` with ``m`` the number of moved
positions and taking ``elem_sym_poly`` over the untouched ones is wrong already at
``p == k``.  The non-equivariant rule ``groth_pieri_mul`` grades instead by
``beta^{d - (number of marked steps)}`` with the multiplicity counting admissible
markings of the chain (``elem_sym_chains_groth``), so the equivariant coefficient
presumably needs that marking data rather than the moved/untouched split.

``var_x`` is unused (kept for signature compatibility).

<a id="schubmult.mult.groth_double.grothmult_double_top"></a>

#### grothmult\_double\_top

```python
def grothmult_double_top(coeff_dict, k, zvar=None, var2=None, beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by the top linear block ``prod_{i=1}^{k}(x_i + zvar)``.

Closed positive Molev--Sagan Pieri rule (conjectural; verified exhaustively on
``S_4`` and sampled through ``S_6``, ``k <= 5``, against
``grothmult_double_block(..., zvar=-zvar, fgl=False)``):

.. math::

    \prod_{i=1}^{k}(x_i + z)\,\mathfrak{G}_u(x; y)
        = \sum_{w} \beta^{\,d - k + |Q|} \Bigl(\prod_{i=1}^{k} f_i\Bigr)\,
          \mathfrak{G}_w(x; y),
    \qquad d = \ell(w) - \ell(u),

where the factor ``f_i`` depends on the fate of the window value ``u(i)``:

* ``u(i) = w(i)`` (the set ``Q``):  ``(z(1 + beta*y_{u(i)}) - y_{u(i)}) / (1 + beta*y_{u(i)})``,
  i.e. ``z (+) (-)y_{u(i)}``, the K-theoretic analogue of ``z - y_{u(i)}``;
* ``u(i)`` stays in the window but moves left:  ``1 - beta*zvar``;
* ``u(i)`` exits the window or moves right within it:  ``1/(1 + beta*y_{u(i)})``.

The sum runs over the marked-chain K-Pieri support (``_top_block_support``).
At ``beta = 0`` this collapses to the ``p = k`` Pieri formula for double
Schubert polynomials [Samuel, Theorem 7.1]:
``S_u(x;y) prod(x_i - z) = sum_{u ->_k w} prod_{i in Q}(y_{u(i)} - z) S_w(x;y)``
with ``z -> -z``.

<a id="schubmult.mult.groth_double.mult_poly_groth_double"></a>

#### mult\_poly\_groth\_double

```python
def mult_poly_groth_double(coeff_dict,
                           poly,
                           var_x=None,
                           var_y=None,
                           beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var_y)`` by an arbitrary polynomial in ``var_x``.

Mirrors ``mult_poly_double``; the leaves of the ``Add``/``Mul``/``Pow`` recursion
are handled by ``single_variable_groth``.

<a id="schubmult.mult.groth_double.dgroth_to_dschub"></a>

#### dgroth\_to\_dschub

```python
def dgroth_to_dschub(v, var3, beta=None)
```

Expand ``G_v(x, var3)`` in double Schubert polynomials: ``{v': coeff}``.

``sum_{v'} coeff_{v'} S_{v'}(x, var3) = G_v(x, var3)`` with coefficients in
``var3`` and ``beta``.  Exact but slow; delegates to ``grothendieck_poly``
with ``keep_as_schub=True``.

<a id="schubmult.mult.groth_double.groth_elem_sym_func"></a>

#### groth\_elem\_sym\_func

```python
def groth_elem_sym_func(k, i, u1, u2, v1, v2, vdiff, varl1, varl2, beta)
```

Expression form of ``_groth_elem_sym_frac``; see there for the rule.

<a id="schubmult.mult.groth_double.grothmult_double"></a>

#### grothmult\_double

```python
def grothmult_double(perm_dict, v, var2=None, var3=None, beta=None)
```

Multiply double Grothendieck polynomials, mirroring ``schubmult_double``.

Computes the expansion of ``sum_u coeff_u G_u(x, var2) * G_v(x, var3)`` in
the basis ``{G_w(x, var2)}`` and returns it as ``{w: coeff_w}``.

``v = s_k`` uses the verified chain formula of Corollary 8.2, and
``max_descent == 1`` folds that column by column.  General ``v`` goes through
``dgroth_to_dschub`` (exact, slow) and the conjectural vpath kernel
``_groth_schub_vpath_mul``, one run per double Schubert ``S_{v'}`` in the
expansion of ``G_v``.

The chain rank is inferred from the current permutation and selected positions.

<a id="schubmult.mult.positivity"></a>

# schubmult.mult.positivity

Manifestly positive representations of double Schubert structure constants.

The structure constants ``c^w_{u,v}(y, z)`` of double Schubert polynomial
multiplication (``schubmult_double``) are known to be polynomials in the
differences ``y_i - z_j`` with nonnegative integer coefficients (Graham's
positivity theorem). This module computes that manifestly positive form:

- ``posify``: the main recursive engine. Reduces ``(u, v, w)`` via known
  combinatorial identities (pattern-avoidance checks, dominance/one-dominance,
  descent and coefficient reductions in ``schubmult.utils.schub_lib``) down to
  cases with closed positive formulas (``dualcoeff``, ``forwardcoeff``, or a
  single elementary symmetric polynomial), falling back to the integer-LP
  solver ``compute_positive_rep`` when no reduction applies.
- ``compute_positive_rep``: expresses an arbitrary such polynomial as a
  nonnegative-integer combination of product-of-differences monomials, found
  via an integer program (PuLP) over a spanning set of candidate monomials.
- ``dualcoeff``/``forwardcoeff``/``dualpieri``: closed-form positive rules for
  special cases (``u`` dominates ``w``, the ``will_formula_work`` forward Monk
  case, and the dual Pieri expansion respectively).

<a id="schubmult.mult.positivity.compute_positive_rep"></a>

#### compute\_positive\_rep

```python
def compute_positive_rep(val, var2=None, var3=None, msg=False)
```

Express ``val`` as a nonnegative-integer combination of product-of-differences monomials.

``val`` must be a polynomial in ``var2``/``var3`` known (by positivity of
double Schubert structure constants) to admit an expansion
``sum_b n_b * prod (var2_i - var3_j)`` with ``n_b >= 0`` integers. Builds a
candidate spanning set of such product monomials from ``val``'s own
monomials, then solves an integer program (via PuLP) for nonnegative
integer coefficients ``n_b`` matching ``val`` exactly.

**Arguments**:

- `val` - Symbolic polynomial expression in ``var2``/``var3``.
- `var2` - First secondary alphabet (``y``).
- `var3` - Second secondary alphabet (``z``).
- `msg` - Passed through to the LP solver as its ``msg`` (verbosity) option.
  

**Returns**:

  A symbolic expression equal to ``val``, written as a sum of
  nonnegative-integer multiples of product-of-differences monomials.
  

**Raises**:

- `Exception` - If the reconstructed expression does not equal ``val``
  (i.e. no valid nonnegative integer solution reproduces it exactly).

<a id="schubmult.mult.positivity.posify"></a>

#### posify

```python
@cached(
    cache={},
    key=lambda val, u2, v2, w2, var2=None, var3=
    None, msg=False, sign_only=False, optimize=True: hashkey(
        val, u2, v2, w2, var2, var3, msg, sign_only, optimize),
)
def posify(val,
           u2,
           v2,
           w2,
           var2=None,
           var3=None,
           msg=False,
           sign_only=False,
           optimize=True,
           n=_vars.n)
```

Manifestly positive representation of the structure constant ``c^{w2}_{u2,v2}(var2, var3)``.

``val`` is the (already computed, possibly not manifestly positive) value of
the coefficient of ``S_{w2}`` in ``S_{u2}(x, var2) * S_{v2}(x, var3)``.
Recursively reduces ``(u2, v2, w2)`` via pattern-avoidance-guarded identities
(``try_reduce_u``/``try_reduce_v``, ``reduce_descents``, ``reduce_coeff``,
``is_split_two``) toward cases handled by closed positive formulas
(a single elementary symmetric polynomial when ``v`` has one nonzero code
entry, ``dualcoeff`` when ``will_formula_work(v, u)`` or ``u`` dominates
``w``, ``forwardcoeff`` when ``will_formula_work(u, v)``, or the
length-one-difference case built from ``pull_out_var``/``schubpoly``
directly). Falls back to ``compute_positive_rep`` (an integer-LP search)
when no reduction or closed formula applies and ``optimize`` is true.

Results are cached by ``(val, u2, v2, w2, var2, var3, msg, sign_only, optimize)``.

**Arguments**:

- `val` - The structure constant to re-express positively.
- `u2` - First factor's permutation.
- `v2` - Second factor's permutation.
- `w2` - Target permutation (coefficient of ``S_{w2}``).
- `var2` - First secondary alphabet.
- `var3` - Second secondary alphabet.
- `msg` - Verbosity flag passed down to ``compute_positive_rep``'s LP solver.
- `sign_only` - If ``True``, only determine and return the sign of ``val``
  (``-1``, ``0``, or ``1``) rather than a full positive expression.
- `optimize` - If ``False``, return ``val`` unchanged when no closed-form
  reduction applies (skip the LP fallback); if ``None`` and that
  case is reached, raise.
- `n` - Size of the ambient alphabet used when no other bound is available.
  

**Returns**:

  A manifestly positive expression equal to ``val`` (or, if
  ``sign_only``, one of ``-1``, ``0``, ``1``).

<a id="schubmult.mult.positivity.shiftsub"></a>

#### shiftsub

```python
def shiftsub(pol, var2=None)
```

Shift every ``var2[i]`` in ``pol`` up to ``var2[i + 1]`` (for ``i`` in ``0..98``).

<a id="schubmult.mult.positivity.posify_generic_partial"></a>

#### posify\_generic\_partial

```python
def posify_generic_partial(val, u2, v2, w2)
```

``posify`` specialized to the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``.

Asserts (raises on mismatch) that the recomputed positive expression equals
the input ``val``, as a consistency check.

<a id="schubmult.mult.positivity.schubmult_generic_partial_posify"></a>

#### schubmult\_generic\_partial\_posify

```python
@cache
def schubmult_generic_partial_posify(u2, v2)
```

Manifestly positive expansion of ``S_{u2}(x, var_g1) * S_{v2}(x, var_g2)``.

Returns ``{w2: coeff}`` where each ``coeff`` is the positive representation
(via ``posify_generic_partial``) of the corresponding
``schubmult_double_pair_generic_alt`` coefficient.

<a id="schubmult.mult.positivity.forwardcoeff"></a>

#### forwardcoeff

```python
def forwardcoeff(u, v, perm, var2=None, var3=None)
```

Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "forward" case
(used when ``will_formula_work(u, v)`` holds in ``posify``).

Writes ``muv = uncode(v.theta())`` and reduces to a lookup in
``schubmult_double_pair(u, muv, var2, var3)`` when the length condition
``(perm * (~v * muv)).inv == (~v * muv).inv + perm.inv`` holds; returns 0 otherwise.

<a id="schubmult.mult.positivity.dualcoeff"></a>

#### dualcoeff

```python
def dualcoeff(u, v, perm, var2=None, var3=None)
```

Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "dual" case
(used in ``posify`` when ``will_formula_work(v, u)`` holds or ``u`` dominates ``perm``).

When ``u`` is the identity, reduces directly to a single Schubert
polynomial ``schubpoly(v * (~perm), var2, var3)``. Otherwise expands via
``dualpieri`` (directly if ``u`` dominates ``perm``, or after rewriting to
``u``'s dominant permutation ``uncode(u.theta())`` otherwise), summing
products of ``(var2_{i+1} - var3_j)`` factors against a final
``schubpoly`` term.

<a id="schubmult.mult.positivity.dualpieri"></a>

#### dualpieri

```python
def dualpieri(mu, v, w)
```

Dual Pieri expansion used by ``dualcoeff``: enumerate the data witnessing
``S_mu * S_v -> S_w`` when ``mu`` is dominant.

Compares ``mu``'s inverse code against ``w``'s inverse code layer by layer,
peeling one "cycle" of variables per layer via ``divdiffable``/``pull_out_var``,
and returns the list of ``[vlist, vp]`` pairs consumed by ``dualcoeff`` to
build the final positive expression (empty list if ``w`` is not reachable
from ``mu``, ``v`` this way).

<a id="schubmult.mult.quantum"></a>

# schubmult.mult.quantum

Quantum (single) Schubert polynomial multiplication.

Implements ``schubmult_q``/``schubmult_q_fast``: the product of a linear
combination of quantum Schubert polynomials ``S_u`` with a single ``S_v``,
returned as a coefficient dict ``{w: coeff}`` polynomial in the quantum
parameters ``q_1, q_2, ...``. Uses the same ``theta``/v-path recursion as
``schubmult.mult.single``, with the elementary-symmetric step generalized to
the quantum Pieri-type moves of ``elem_sym_perms_q`` (which may pick up a
factor of ``q`` when a Bruhat move is replaced by its quantum analogue).
``schubmult_q`` uses ``strict_theta`` (no repeated layer merging);
``schubmult_q_fast``/``_schubmult_q_fast_python`` uses ``medium_theta`` and
merges adjacent equal-length layers via ``double_elem_sym_q`` for speed.

<a id="schubmult.mult.quantum.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var_q=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum`` (quantum Monk rule).

Same structure as ``schubmult.mult.single.single_variable``, using
``elem_sym_perms_q`` so that some Bruhat moves carry a factor from ``var_q``.

<a id="schubmult.mult.quantum.mult_poly_q"></a>

#### mult\_poly\_q

```python
def mult_poly_q(coeff_dict, poly, var_x=_vars.var_x, var_q=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_py``.

<a id="schubmult.mult.quantum.schubmult_q_fast"></a>

#### schubmult\_q\_fast

```python
def schubmult_q_fast(perm_dict, v, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available, falling
back to ``_schubmult_q_fast_python`` (the ``medium_theta``-based recursion
with merged equal-length layers) otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the quantum Schubert
  polynomial to multiply by.
- `q_var` - Generating set for the quantum parameters.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``, polynomial in ``q_var``.

<a id="schubmult.mult.quantum.schubmult_q"></a>

#### schubmult\_q

```python
def schubmult_q(perm_dict, v)
```

Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

Reference (non-"fast") implementation: uses ``strict_theta`` and processes
every layer individually (no merging of equal-length adjacent layers), so it
is simpler but slower than ``schubmult_q_fast``. Results agree with
``schubmult_q_fast`` for all inputs.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the quantum Schubert
  polynomial to multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``, polynomial in the
  default quantum parameters ``q``.

<a id="schubmult.mult.quantum_double"></a>

# schubmult.mult.quantum\_double

Quantum double Schubert polynomial multiplication.

Implements ``schubmult_q_double``/``schubmult_q_double_fast``: the product of a
linear combination of quantum double Schubert polynomials ``S_u(x, var2)`` with
a single ``S_v(x, var3)``, returned as a coefficient dict ``{w: coeff}``
polynomial in ``var2``, ``var3``, and the quantum parameters ``q_1, q_2, ...``.
Uses the same ``theta``/v-path recursion as ``schubmult.mult.double``, with the
elementary-symmetric step generalized to the quantum moves of
``elem_sym_perms_q`` and the coefficient function to ``elem_sym_func_q``.

Also provides: ``mult_poly_q_double`` (multiply by an arbitrary polynomial),
``apply_peterson_woodward`` (parabolic quantum via the Peterson-Woodward
comparison theorem), ``q_posify``/``q_partial_posify_generic`` (manifestly
positive display of quantum structure constants), ``schubpoly_quantum``
(the quantum Schubert polynomial itself), ``nil_hecke`` (quantum nilHecke
action), and ``factor_out_q`` (split a polynomial by its ``q``-monomials).

<a id="schubmult.mult.quantum_double.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var_y=None, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by the single variable ``x_varnum`` (quantum equivariant Monk rule).

The diagonal term contributes ``var_y[u(varnum)]``; the off-diagonal terms come from
``elem_sym_positional_perms_q``, each carrying its ``q``-monomial and sign.

<a id="schubmult.mult.quantum_double.mult_poly_q_double"></a>

#### mult\_poly\_q\_double

```python
def mult_poly_q_double(coeff_dict,
                       poly,
                       var_x=None,
                       var_y=None,
                       q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_double``.

<a id="schubmult.mult.quantum_double.mult_poly_q_double_alt"></a>

#### mult\_poly\_q\_double\_alt

```python
def mult_poly_q_double_alt(coeff_dict,
                           poly,
                           var_x=None,
                           var_y=None,
                           q_var=_vars.q_var)
```

Variant of ``mult_poly_q_double`` that folds each factor via ``schubmult_q_double_dict_fast``
instead of ``single_variable``.

<a id="schubmult.mult.quantum_double.nil_hecke"></a>

#### nil\_hecke

```python
def nil_hecke(perm_dict, v, n, var2=None, var3=None)
```

Quantum nilHecke action: like ``schubmult_q_double`` but using the descent-side
``elem_sym_perms_q_op`` moves (bounded by ``n``) with ``up`` and ``up2`` swapped in the
coefficient function.

<a id="schubmult.mult.quantum_double.schubmult_q_double_pair"></a>

#### schubmult\_q\_double\_pair

```python
@cache
def schubmult_q_double_pair(perm1, perm2, var2=None, var3=None, q_var=None)
```

``schubmult_q_double_fast`` specialized to a single ``perm1`` with coefficient 1, cached.

<a id="schubmult.mult.quantum_double.schubmult_q_double_pair_generic"></a>

#### schubmult\_q\_double\_pair\_generic

```python
@cache
def schubmult_q_double_pair_generic(perm1, perm2)
```

``schubmult_q_double_pair`` with the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``/``_vars.q_var``.

<a id="schubmult.mult.quantum_double.schubmult_q_generic_partial_posify"></a>

#### schubmult\_q\_generic\_partial\_posify

```python
@cache
def schubmult_q_generic_partial_posify(u2, v2)
```

Manifestly positive (where possible) expansion of ``S_{u2} * S_{v2}`` over the generic alphabets,
applying ``q_partial_posify_generic`` to each coefficient.

<a id="schubmult.mult.quantum_double.q_posify"></a>

#### q\_posify

```python
def q_posify(u, v, w, val, var2, var3, q_var, msg)
```

Manifestly positive representation of the quantum double structure constant ``c^w_{u,v}``.

Splits ``val`` by ``q``-monomial (``factor_out_q``), then for each piece either takes it
as-is (integer, or when ``v``'s inverse code is already in medium-theta form), reduces the
triple ``(u, v, w)`` via ``reduce_q_coeff`` until the ``q``-monomial becomes trivial and
delegates to the classical ``posify``, or falls back to ``compute_positive_rep``.
Raises if the reconstruction does not equal ``val``.

<a id="schubmult.mult.quantum_double.q_partial_posify_generic"></a>

#### q\_partial\_posify\_generic

```python
def q_partial_posify_generic(val, u, v, w)
```

Like ``q_posify`` over the generic alphabets, but only attempts positivity when ``v`` contains
a ``1432`` or ``312`` pattern (otherwise the raw value is already manifestly positive), and
leaves non-reducible ``q``-pieces unchanged rather than running the LP.

<a id="schubmult.mult.quantum_double.apply_peterson_woodward"></a>

#### apply\_peterson\_woodward

```python
def apply_peterson_woodward(coeff_dict, parabolic_index, q_var=_vars.q_var)
```

Project a full-flag quantum product onto the parabolic quantum cohomology for ``parabolic_index``.

Implements the Peterson-Woodward comparison: for each ``q``-monomial of each coefficient,
checks the ``omega``/``check_blocks`` compatibility conditions on the exponent vector,
multiplies the indexing permutation by the appropriate parabolic longest elements, keeps
only the ``parabolic``-minimal results, and reindexes the surviving ``q`` variables.

**Arguments**:

- `coeff_dict` - Full-flag quantum coefficient dict ``{Permutation: coeff}``.
- `parabolic_index` - Sorted list of 1-indexed positions generating the parabolic subgroup.
- `q_var` - Quantum parameter generating set.
  

**Returns**:

- `dict` - Parabolic quantum coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.quantum_double.elem_sym_func_q_q"></a>

#### elem\_sym\_func\_q\_q

```python
def elem_sym_func_q_q(k,
                      i,
                      u1,
                      u2,
                      v1,
                      v2,
                      udiff,
                      vdiff,
                      varl1,
                      varl2,
                      q_var=_vars.q_var)
```

Fully-quantum coefficient function for the v-path recursion (used by ``schubpoly_quantum``):
the quantum elementary symmetric polynomial ``elem_sym_poly_q`` in the fixed-window ``y``
variables and the ``call_zvars`` ``z`` variables.

<a id="schubmult.mult.quantum_double.schubpoly_quantum"></a>

#### schubpoly\_quantum

```python
def schubpoly_quantum(v, var_x=None, var_y=None, q_var=_vars.q_var, coeff=1)
```

The quantum double Schubert polynomial ``S_v(var_x, var_y)`` itself, as a symbolic expression.

Runs the v-path recursion starting from the identity with ``elem_sym_func_q_q`` and reads off
the coefficient of the identity permutation.

<a id="schubmult.mult.quantum_double.schubmult_q_double"></a>

#### schubmult\_q\_double

```python
def schubmult_q_double(perm_dict, v, var2=None, var3=None, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the quantum double Schubert polynomial ``S_v(x, var3)``.

Reference (non-"fast") implementation: uses ``strict_theta`` and processes every layer
individually. Results agree with ``schubmult_q_double_fast``.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
- `q_var` - Quantum parameter generating set.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.quantum_double.schubmult_q_double_dict_fast"></a>

#### schubmult\_q\_double\_dict\_fast

```python
def schubmult_q_double_dict_fast(perm_dict1,
                                 perm_dict2,
                                 var2=None,
                                 var3=None,
                                 q_var=_vars.q_var)
```

Multiply two coefficient dicts of quantum double Schubert polynomials together.

Sums ``schubmult_q_double_fast(perm_dict1, v, ...)`` scaled by ``coeff2_v`` over ``v`` in ``perm_dict2``.

<a id="schubmult.mult.quantum_double.schubmult_q_double_fast"></a>

#### schubmult\_q\_double\_fast

```python
def schubmult_q_double_fast(perm_dict,
                            v,
                            var2=None,
                            var3=None,
                            q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the quantum double Schubert polynomial ``S_v(x, var3)``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available (and both secondary
alphabets are given), falling back to ``_schubmult_q_double_fast_python`` (the
``medium_theta``-based recursion with merged equal-length layers) otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
- `q_var` - Quantum parameter generating set.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.quantum_double.sum_q_dict"></a>

#### sum\_q\_dict

```python
def sum_q_dict(q_dict1, q_dict2)
```

Add two ``{q_monomial: coeff}`` dicts.

<a id="schubmult.mult.quantum_double.mul_q_dict"></a>

#### mul\_q\_dict

```python
def mul_q_dict(q_dict1, q_dict2)
```

Multiply two ``{q_monomial: coeff}`` dicts (convolution over monomials).

<a id="schubmult.mult.quantum_double.factor_out_q"></a>

#### factor\_out\_q

```python
def factor_out_q(poly, q_var=_vars.q_var)
```

Split ``poly`` by its ``q``-monomials: return ``{q_monomial: coefficient}`` with coefficients
free of ``q_var`` variables. Recurses over the ``Add``/``Mul``/``Pow`` structure; a polynomial
with no ``q`` variables maps to ``{1: poly}``.

<a id="schubmult.mult.separated_descents"></a>

# schubmult.mult.separated\_descents

Separated-descents product of (double) Grothendieck polynomials.

Implements the *pipe puzzle* formula of Fan--Guo--Xiong,
"Bumpless pipe dreams meet puzzles" (arXiv:2309.00467), Theorem 2.5.

For permutations ``u`` and ``v`` with *separated descents* at a position ``k``
(i.e. every descent of ``u`` is ``<= k`` and every descent of ``v`` is ``>= k``)
the double Grothendieck polynomials satisfy

    G_u(x, y) * G_v(x, t) = sum_w c_{u,v}^w(t, y) * G_w(x, t),

and the structure constants ``c_{u,v}^w(t, y)`` are given by a positive
(in the sense of Theorem 2.5) sum over pipe puzzles.

The beta convention matches the rest of ``schubmult``: the multiplicative
formal group law is ``x (+) y = x + y + beta*x*y`` with formal subtraction

    x (-) y = (x - y) / (1 + beta*y).

Setting ``beta = 0`` recovers the double Schubert (cohomology) structure
constants of Theorem 4.x (the "Schubert pipe puzzle" specialization).

<a id="schubmult.mult.separated_descents.separated_descents_coeffs"></a>

#### separated\_descents\_coeffs

```python
def separated_descents_coeffs(u, v, var1, var2, beta=None, grid_size=None)
```

Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
``var2`` are the ``y`` variables (secondary variables of ``u``).

Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

In K-theory the product may involve ``G_w`` with ``w`` in a larger
symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
``grid_size`` is not supplied it is chosen large enough to capture every
such ``w``: the maximum value of any appearing ``w`` is bounded by
``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.

<a id="schubmult.mult.separated_descents.separated_descents_coeffs_plus"></a>

#### separated\_descents\_coeffs\_plus

```python
def separated_descents_coeffs_plus(u,
                                   v,
                                   var1,
                                   var2,
                                   beta=None,
                                   grid_size=None,
                                   mangle_genset=False)
```

Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
``var2`` are the ``y`` variables (secondary variables of ``u``).

Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

In K-theory the product may involve ``G_w`` with ``w`` in a larger
symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
``grid_size`` is not supplied it is chosen large enough to capture every
such ``w``: the maximum value of any appearing ``w`` is bounded by
``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.

<a id="schubmult.mult.separated_descents.grothmult_double"></a>

#### grothmult\_double

```python
def grothmult_double(perm_dict, v, var1, var2, beta=None)
```

Separated-descents product of double Grothendieck polynomials.

Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

    coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.

<a id="schubmult.mult.separated_descents.grothmult_double_plus"></a>

#### grothmult\_double\_plus

```python
def grothmult_double_plus(perm_dict,
                          v,
                          var1,
                          var2,
                          beta=None,
                          mangle_genset=False)
```

Separated-descents product of double Grothendieck polynomials.

Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

    coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.

<a id="schubmult.mult.single"></a>

# schubmult.mult.single

Ordinary (single) Schubert polynomial multiplication.

Implements the ``schubmult_py`` kernel: the product of a linear combination of
Schubert polynomials ``S_u`` (given as a dict ``{u: coeff}``) with a single
Schubert polynomial ``S_v``, returned as a coefficient dict ``{w: coeff}``.

The algorithm is the recursive "v-path" / transition method used throughout
``schubmult``: write ``theta = (~v).theta()`` for the dominant weakly-decreasing
vector bounding ``v``'s Lehmer code, let ``mu = uncode(theta)`` and
``vmu = v * mu``, and process the entries of ``theta`` one at a time, tracking
Bruhat-chain "v-paths" from ``vmu`` down to the identity (``compute_vpathdicts``)
alongside chains of elementary-symmetric moves on the ``u`` side
(``elem_sym_perms``). Accumulating consistent pairs of chains and reading off
the coefficient landing on ``vmu`` gives the product.

<a id="schubmult.mult.single.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum)
```

Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum``.

Uses the classical Monk rule: ``x_k * S_u = sum S_{u t_{ij}}`` over Bruhat
covers ``u t_{ij}`` with ``i <= k < j`` (added) minus those with ``j <= k < i``
(subtracted), via ``elem_sym_perms(u, 1, varnum)``.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `varnum` - 1-indexed variable index ``k``.
  

**Returns**:

- `dict` - The updated coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.mult_poly_py"></a>

#### mult\_poly\_py

```python
def mult_poly_py(coeff_dict, poly, var_x=_vars.var_x)
```

Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``; each single
variable leaf is dispatched to ``single_variable``, and any other leaf just
scales every coefficient.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `poly` - Symbolic polynomial expression in the variables of ``var_x``.
- `var_x` - Generating set identifying the ``x`` variables (default ``x``).
  

**Returns**:

- `dict` - The updated coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.schubmult_py"></a>

#### schubmult\_py

```python
def schubmult_py(perm_dict, v)
```

Multiply ``sum_u coeff_u S_u(x)`` by the (ordinary) Schubert polynomial ``S_v``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available and the
permutations fit within its ``MAXN``, falling back to the pure-Python
implementation otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}`` (integer coefficients).
- `v` - Permutation (or array-form list) indexing the Schubert polynomial to
  multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}`` for the product.

<a id="schubmult.mult.single.schubmult_py_down"></a>

#### schubmult\_py\_down

```python
def schubmult_py_down(perm_dict, v)
```

Divided-difference ("down") variant of ``_schubmult_py_python``.

Same v-path recursion but built from ``elem_sym_perms_op`` (Bruhat *descents*)
instead of ``elem_sym_perms``, used for the down/dual side of the transition
recursion rather than ordinary multiplication.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.schub_coprod_py"></a>

#### schub\_coprod\_py

```python
def schub_coprod_py(perm, indices)
```

Coproduct of ``S_perm`` restricted to the variable split named by ``indices``.

Computes the expansion of the (single) Schubert polynomial coproduct
``Delta_{indices}(S_perm) = sum (firstperm, secondperm) -> coeff`` by
multiplying the Grassmannian permutation for ``indices`` against ``perm``
via ``schubmult_py`` and splitting each resulting permutation's window into
its first ``N`` and remaining ``len(perm) - N`` values.

**Arguments**:

- `perm` - Permutation (or array-form list) to take the coproduct of.
- `indices` - Iterable of 1-indexed positions selecting the variable split.
  

**Returns**:

- `dict` - Mapping ``{(firstperm, secondperm): coeff}``.

<a id="schubmult.rings"></a>

# schubmult.rings

Ring structures built on combinatorial bases.

Subpackages:

- `schubert`: the Schubert-family rings (``Sx``, ``DSx``, ``Gx``, ``QSx``, ...) -- the
  main user-facing algebra interface.
- `polynomial_algebra`: the polynomial ring ``Z[x_1, x_2, ...]`` with pluggable bases
  (monomial, Schubert, key, slide, Grothendieck, ...).
- `free_algebra`: noncommutative free algebras with combinatorial bases.
- `combinatorial`: rings whose basis elements are combinatorial objects (RC graphs,
  BPDs, plactic classes, ...) rather than permutations.

Modules here: `base_ring` (the shared dict-based ring/element machinery), `tensor_ring`,
`product_ring`, `direct_product_ring`, `printing` (display symbols), `nsym`,
`quasisymmetric_functions`, `thompson_algebra`.

This ``__init__`` re-exports nothing (the commented-out block below is legacy);
import from the subpackages directly.

<a id="schubmult.rings.base_ring"></a>

# schubmult.rings.base\_ring

Shared machinery for every ring in `schubmult.rings`.

`BaseRingElement` is a ``dict`` mapping basis keys (permutations, RC graphs, tuples,
...) to coefficients, wired into sympy's printing and arithmetic protocols so that
elements can be added, multiplied, and displayed. `BaseRing` provides the
corresponding ring-level operations (``add``/``sub``/``mul``, coercion via
``domain_new``, construction via ``from_dict``/``from_expr``) and declares the hooks
concrete rings must implement (``new``, ``printing_term``, ``mul_expr``, ...).
A ring's element type is created dynamically as ``self.dtype`` with ``ring`` bound.

<a id="schubmult.rings.base_ring.BaseRingElement"></a>

## BaseRingElement Objects

```python
class BaseRingElement(DomainElement, DefaultPrinting, dict)
```

A ring element: ``{basis_key: coefficient}`` with sympy-compatible arithmetic and printing.

<a id="schubmult.rings.base_ring.BaseRingElement.is_zero"></a>

#### is\_zero

```python
@property
def is_zero()
```

Whether every coefficient is exactly zero.

<a id="schubmult.rings.base_ring.BaseRingElement.parent"></a>

#### parent

```python
def parent()
```

The ring this element belongs to (sympy domain protocol).

<a id="schubmult.rings.base_ring.BaseRingElement.has_free"></a>

#### has\_free

```python
def has_free(*args)
```

Whether any of the given symbols appears in ``free_symbols``.

<a id="schubmult.rings.base_ring.BaseRingElement.apply_to_keys"></a>

#### apply\_to\_keys

```python
def apply_to_keys(func)
```

Map each basis key through ``func`` (dropping keys where it returns ``None``), keeping coefficients.

<a id="schubmult.rings.base_ring.BaseRingElement.as_terms"></a>

#### as\_terms

```python
def as_terms()
```

Terms ``coeff * basis_symbol`` in dict order (sympy printing hook).

<a id="schubmult.rings.base_ring.BaseRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by basis key (sympy printing hook).

<a id="schubmult.rings.base_ring.BaseRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Coproduct into the tensor square ring, via ``ring.coproduct_on_basis``.

<a id="schubmult.rings.base_ring.BaseRingElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

``{basis_symbol: coeff}`` mapping display symbols to coefficients.

<a id="schubmult.rings.base_ring.BaseRingElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

With ``deep=True`` expand to an explicit polynomial (``as_polynomial``); with ``deep=False`` only
expand each coefficient, keeping the basis.

<a id="schubmult.rings.base_ring.BaseRingElement.as_expr"></a>

#### as\_expr

```python
def as_expr()
```

Sum of the ``as_terms()`` as a sympy ``Add``.

<a id="schubmult.rings.base_ring.BaseRingElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Hook: expand this element to an explicit polynomial expression.

<a id="schubmult.rings.base_ring.BaseRingElement.almosteq"></a>

#### almosteq

```python
def almosteq(other)
```

Equality up to coefficient expansion.

<a id="schubmult.rings.base_ring.BaseRingElement.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

Tensor product ``self (x) other`` in the `TensorRing` of the two rings.

<a id="schubmult.rings.base_ring.BaseRing"></a>

## BaseRing Objects

```python
class BaseRing(Ring, CompositeDomain)
```

Abstract base ring over a sympy coefficient domain (default ``EXRAW``); see the module docstring.

<a id="schubmult.rings.base_ring.BaseRing.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

The `TensorRing` ``self (x) other``.

<a id="schubmult.rings.base_ring.BaseRing.to_sympy"></a>

#### to\_sympy

```python
def to_sympy(elem)
```

Convert an element to a sympy expression (``as_expr``).

<a id="schubmult.rings.base_ring.BaseRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(domain=None)
```

**Arguments**:

- `domain` - Coefficient domain; defaults to sympy's ``EXRAW`` (arbitrary expressions).

<a id="schubmult.rings.base_ring.BaseRing.add"></a>

#### add

```python
def add(elem, other)
```

Coefficient-wise sum, dropping zeros.

<a id="schubmult.rings.base_ring.BaseRing.sub"></a>

#### sub

```python
def sub(elem, other)
```

Coefficient-wise difference, dropping zeros.

<a id="schubmult.rings.base_ring.BaseRing.neg"></a>

#### neg

```python
def neg(elem)
```

Negate every coefficient.

<a id="schubmult.rings.base_ring.BaseRing.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Right-multiply by a scalar (via ``domain_new``), falling back to ``mul_expr``.

<a id="schubmult.rings.base_ring.BaseRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply by a scalar (via ``domain_new``), falling back to ``mul_expr``.

<a id="schubmult.rings.base_ring.BaseRing.from_sympy"></a>

#### from\_sympy

```python
def from_sympy(expr)
```

Alias for ``from_expr`` (sympy domain protocol).

<a id="schubmult.rings.base_ring.BaseRing.new"></a>

#### new

```python
def new(x)
```

Hook: build an element from ``x``.

<a id="schubmult.rings.base_ring.BaseRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Hook: the sympy symbol displayed for basis key ``k``.

<a id="schubmult.rings.base_ring.BaseRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Hook: coproduct of basis key ``k`` in the tensor square ring.

<a id="schubmult.rings.base_ring.BaseRing.one"></a>

#### one

```python
@property
def one()
```

The multiplicative identity: coefficient 1 on ``zero_monom``.

<a id="schubmult.rings.base_ring.BaseRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(elem)
```

Hook: whether ``elem`` should use the ``elem_mul`` fast path.

<a id="schubmult.rings.base_ring.BaseRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Hook: elementary-symmetric fast-path multiplication.

<a id="schubmult.rings.base_ring.BaseRing.from_dict"></a>

#### from\_dict

```python
def from_dict(element, orig_domain=None)
```

Build an element from ``{key: coeff}``, coercing each coefficient via ``domain_new`` and dropping zeros.

<a id="schubmult.rings.base_ring.BaseRing.from_dict_unchecked"></a>

#### from\_dict\_unchecked

```python
def from_dict_unchecked(element)
```

from_dict for coefficients already known to lie in the domain (drops structural zeros only).

<a id="schubmult.rings.base_ring.BaseRing.zero"></a>

#### zero

```python
@property
def zero()
```

The empty element.

<a id="schubmult.rings.base_ring.BaseRing.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce ``element`` into the coefficient domain (``sympify``), refusing ring/domain elements.

<a id="schubmult.rings.base_ring.BaseRing.from_expr"></a>

#### from\_expr

```python
def from_expr(x)
```

Build an element from an expression by multiplying the identity by it.

<a id="schubmult.rings.base_ring.BaseRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Hook: multiply ``elem`` by a symbolic expression ``x``.

<a id="schubmult.rings.combinatorial"></a>

# schubmult.rings.combinatorial

Combinatorial rings: rings whose basis elements are combinatorial objects (RC graphs, BPDs, WC graphs,
tableaux, ...) rather than permutations.

The central object is `RCGraphRing`; most other rings here are quotients or variants of it that snap
products to canonical representatives (`HWRCGraphRing`, `KeyRCGraphRing`, `QYRCGraphRing`,
`ForestRCGraphRing`, `SlideRCGraphRing`, ...). `BoundedRCFactorAlgebra` and `GrassTensorAlgebra`
provide factorizations into Grassmannian pieces used to compute products and coproducts.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring"></a>

# schubmult.rings.combinatorial.alt\_rc\_graph\_ring

`AltRCGraphRing`: an alternate `RCGraphRing` implementation exploring a different polynomial
product (``%``) construction; the exported `RCGraphRing` is the primary one.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement"></a>

## AltRCGraphRingElement Objects

```python
class AltRCGraphRingElement(CrystalGraphRingElement,
                            SchubertMonomialRingElement)
```

AltRCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product (only defined when the right side is a dominant RC graph);
the product * is the dual product, defined for any pair of RC graphs.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.vertical_coproduct"></a>

#### vertical\_coproduct

```python
def vertical_coproduct()
```

Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
and induces the mul product.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an AltRCGraphRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an AltRCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

Linear extension of RCGraph.crystal_reflection:
For each basis RCGraph, apply its crystal_reflection(index) and collect results.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRing"></a>

## AltRCGraphRing Objects

```python
class AltRCGraphRing(SchubertMonomialRing, CrystalGraphRing)
```

Alternate RC graph ring; see the module docstring.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the AltRCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_rc\_factor\_algebra

`BoundedRCFactorAlgebra`: a tensor-like algebra on tuples of full Grassmannian RC graphs of
bounded size, factoring Schubert classes into elementary-symmetric pieces (the CEM basis).

Used as the engine behind `RCGraphRing.coproduct_on_basis`, `SchubertRCGraphRing`, and
`SlideRCGraphRing`: ``schub_elem(perm, length)`` gives the factorization of ``S_perm``,
``key_to_rc_graph`` / ``to_rc_graph_ring_element`` squash a key's factors back into an RC graph.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorPrintingTerm"></a>

## BoundedRCFactorPrintingTerm Objects

```python
class BoundedRCFactorPrintingTerm(PrintingTerm)
```

Display symbol for a `BoundedRCFactorAlgebra` basis key.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebraElement"></a>

## BoundedRCFactorAlgebraElement Objects

```python
class BoundedRCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebraElement.prune"></a>

#### prune

```python
def prune()
```

Merge terms whose evaluated RC graphs share forest/omega invariants.

Terms are grouped by
``(rc.forest_weight, rc.omega_invariant[1])`` where
``rc = self.ring.key_to_rc_graph(key)``.
One key per group is kept (first encountered), and coefficients are summed.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebra"></a>

## BoundedRCFactorAlgebra Objects

```python
class BoundedRCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_rc\_forest\_factor\_algebra

`BoundedRCForestFactorAlgebra`: variant of `BoundedRCFactorAlgebra` whose factorizations are
snapped to forest-class representatives (``_to_forest``), for the forest-polynomial setting.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCFactorPrintingTerm"></a>

## BoundedRCFactorPrintingTerm Objects

```python
class BoundedRCFactorPrintingTerm(PrintingTerm)
```

Display symbol for a `BoundedRCForestFactorAlgebra` basis key.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebraElement"></a>

## BoundedRCForestFactorAlgebraElement Objects

```python
class BoundedRCForestFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCForestFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra"></a>

## BoundedRCForestFactorAlgebra Objects

```python
class BoundedRCForestFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_wc\_factor\_algebra

`BoundedWCFactorAlgebra`: the `WCGraph` (K-theoretic) analogue of `BoundedRCFactorAlgebra`,
factoring Grothendieck classes into tuples of full Grassmannian WC graphs.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorPrintingTerm"></a>

## BoundedWCFactorPrintingTerm Objects

```python
class BoundedWCFactorPrintingTerm(PrintingTerm)
```

Display symbol for a `BoundedWCFactorAlgebra` basis key.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebraElement"></a>

## BoundedWCFactorAlgebraElement Objects

```python
class BoundedWCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedWCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebra"></a>

## BoundedWCFactorAlgebra Objects

```python
class BoundedWCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian WC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
WC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebra.key_to_wc_graph"></a>

#### key\_to\_wc\_graph

```python
def key_to_wc_graph(key) -> WCGraph
```

Evaluate a tensor key to an WCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.bpd_ring"></a>

# schubmult.rings.combinatorial.bpd\_ring

`BPDRing`: a `SchubertMonomialRing` whose basis elements are bumpless pipe dreams (`BPD`),
with conversion to `RCGraphRing` via ``to_rc_graph_ring_element``.

<a id="schubmult.rings.combinatorial.bpd_ring.BPDRingElement"></a>

## BPDRingElement Objects

```python
class BPDRingElement(SchubertMonomialRingElement)
```

Linear combination of `BPD` basis elements.

<a id="schubmult.rings.combinatorial.bpd_ring.BPDRingElement.to_rc_graph_ring_element"></a>

#### to\_rc\_graph\_ring\_element

```python
def to_rc_graph_ring_element(
        rc_ring: RCGraphRing | None = None) -> RCGraphRingElement
```

Convert each BPD to its RC graph and re-express in an `RCGraphRing`.

<a id="schubmult.rings.combinatorial.bpd_ring.BPDRing"></a>

## BPDRing Objects

```python
class BPDRing(SchubertMonomialRing)
```

The ring of bumpless pipe dreams; products use `BPD.product`.

<a id="schubmult.rings.combinatorial.chute_move_ring"></a>

# schubmult.rings.combinatorial.chute\_move\_ring

`ChuteMoveRing`: a `SchubertMonomialRing` whose basis elements are `ChuteMoveElement`s
(RC graphs marked with a set of simultaneous chute-move rows).

<a id="schubmult.rings.combinatorial.chute_move_ring.ChuteMoveRingElement"></a>

## ChuteMoveRingElement Objects

```python
class ChuteMoveRingElement(SchubertMonomialRingElement)
```

ChuteMoveRing elements are linear combinations of ChuteMoveElement basis elements.

<a id="schubmult.rings.combinatorial.chute_move_ring.ChuteMoveRing"></a>

## ChuteMoveRing Objects

```python
class ChuteMoveRing(SchubertMonomialRing)
```

The ring of `ChuteMoveElement`s; products use `ChuteMoveElement.product`.

<a id="schubmult.rings.combinatorial.crystal_graph_ring"></a>

# schubmult.rings.combinatorial.crystal\_graph\_ring

`CrystalGraphRing`: a ring whose basis elements are crystal-graph objects, with the crystal
operators extended linearly to ring elements. Base class for the RC graph / WC graph /
factor-algebra rings in this package.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRing"></a>

## CrystalGraphRing Objects

```python
class CrystalGraphRing(BaseRing)
```

Ring whose basis elements are CrystalGraph-like objects.

We deliberately do not special-case tensor objects here: CrystalGraphTensor
implements the same CrystalGraph API and will be handled by polymorphism.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRing.dtype"></a>

#### dtype

```python
def dtype()
```

A fresh empty element bound to this ring.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement"></a>

## CrystalGraphRingElement Objects

```python
class CrystalGraphRingElement(BaseRingElement, CrystalGraph)
```

Element of the CrystalGraphRing.

Keys are arbitrary objects that implement the CrystalGraph API (including
CrystalGraphTensor). All crystal operators / statistics are lifted linearly
by delegating to the underlying key's methods.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.phi"></a>

#### phi

```python
def phi(index: int) -> int
```

Maximum of ``phi(index)`` over the basis keys.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index: int) -> int
```

Maximum of ``epsilon(index)`` over the basis keys.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index: int)
```

Linearized raising operator: delegate to each key's raising_operator
and collect results in the ring.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index: int)
```

Linearized lowering operator: delegate to each key's lowering_operator.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length() -> int
```

Maximum of ``crystal_length()`` over the basis keys.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight() -> Tuple["CrystalGraphRingElement", Tuple[int, ...]]
```

Apply linearized raising operators until none changes the element; returns ``(element, raise_seq)``.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring"></a>

# schubmult.rings.combinatorial.dual\_rc\_graph\_ring

`DualRCGraphRing`: RC graph ring carrying the dual (polynomial-side) product, computed by
expanding into the Schubert polynomial basis of `PolynomialAlgebra` and back.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement"></a>

## DualRCGraphRingElement Objects

```python
class DualRCGraphRingElement(SchubertMonomialRingElement)
```

DualRCGraphRing elements are linear combinations of RCGraph basis elements under the dual product.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRing"></a>

## DualRCGraphRing Objects

```python
class DualRCGraphRing(SchubertMonomialRing)
```

The dual RC graph ring; see the module docstring.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the DualRCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.eg_plactic_ring"></a>

# schubmult.rings.combinatorial.eg\_plactic\_ring

`EGPlacticRing`: a `CrystalGraphRing` on pairs ``((NilPlactic, length), Plactic)`` -- an RC
graph's Edelman-Greene insertion tableau together with its plactic recording tableau -- with
conversion to and from `RCGraphRing`.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticPrintingTerm"></a>

## EGPlacticPrintingTerm Objects

```python
class EGPlacticPrintingTerm(PrintingTerm)
```

Display symbol for an `EGPlacticRing` basis key.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement"></a>

## EGPlacticRingElement Objects

```python
class EGPlacticRingElement(CrystalGraphRingElement)
```

EGPlacticRing elements are linear combinations of ``((NilPlactic, length), Plactic)`` basis keys.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an EGPlacticRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an EGPlacticRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an EGPlacticRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRing"></a>

## EGPlacticRing Objects

```python
class EGPlacticRing(CrystalGraphRing)
```

The EG-plactic ring; see the module docstring.

<a id="schubmult.rings.combinatorial.eg_ring"></a>

# schubmult.rings.combinatorial.eg\_ring

`EGRing`: a ring whose basis elements are ``(NilPlactic, length)`` pairs -- the Edelman-Greene
insertion tableau of an RC graph's word together with its row count. ``from_rc_graph`` maps an RC
graph to its EG class.

<a id="schubmult.rings.combinatorial.eg_ring.EGPrintingTerm"></a>

## EGPrintingTerm Objects

```python
class EGPrintingTerm(PrintingTerm)
```

Display symbol for an `EGRing` basis key (prints the key directly).

<a id="schubmult.rings.combinatorial.eg_ring.EGRingElement"></a>

## EGRingElement Objects

```python
class EGRingElement(BaseRingElement)
```

Linear combination of ``(NilPlactic, length)`` basis keys.

<a id="schubmult.rings.combinatorial.eg_ring.EGRing"></a>

## EGRing Objects

```python
class EGRing(BaseRing)
```

The Edelman-Greene tableau ring; see the module docstring.

<a id="schubmult.rings.combinatorial.forest_invariant_rc_ring"></a>

# schubmult.rings.combinatorial.forest\_invariant\_rc\_ring

`ForestInvariantRCGraphRing`: `RCGraphRing` quotient where every product is snapped to its crystal
highest-weight representative (identical in behavior to `HWRCGraphRing`).

<a id="schubmult.rings.combinatorial.forest_invariant_rc_ring.ForestInvariantRCGraphRing"></a>

## ForestInvariantRCGraphRing Objects

```python
class ForestInvariantRCGraphRing(RCGraphRing)
```

`RCGraphRing` with products projected onto highest-weight RC graphs (``_snap_highest_weight``).

<a id="schubmult.rings.combinatorial.forest_rc_ring"></a>

# schubmult.rings.combinatorial.forest\_rc\_ring

`ForestRCGraphRing` / `DualForestRCGraphRing`: `RCGraphRing` quotients modeling forest polynomials
(Nadeau-Spink-Tewari). RC graphs are snapped to canonical representatives of their
``forest_weight`` class; the dual variant carries the dual product.

<a id="schubmult.rings.combinatorial.forest_rc_ring.ForestRCGraphRingElement"></a>

## ForestRCGraphRingElement Objects

```python
class ForestRCGraphRingElement(RCGraphRingElement)
```

Element of `ForestRCGraphRing`.

<a id="schubmult.rings.combinatorial.forest_rc_ring.DualForestRCGraphRingElement"></a>

## DualForestRCGraphRingElement Objects

```python
class DualForestRCGraphRingElement(ForestRCGraphRingElement)
```

Element of `DualForestRCGraphRing`.

<a id="schubmult.rings.combinatorial.forest_rc_ring.ForestRCGraphRing"></a>

## ForestRCGraphRing Objects

```python
class ForestRCGraphRing(RCGraphRing)
```

`RCGraphRing` snapped to forest-class representatives; see the module docstring.

<a id="schubmult.rings.combinatorial.forest_rc_ring.DualForestRCGraphRing"></a>

## DualForestRCGraphRing Objects

```python
class DualForestRCGraphRing(RCGraphRing)
```

Dual-product variant of `ForestRCGraphRing`.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra"></a>

# schubmult.rings.combinatorial.grass\_tensor\_algebra

`GrassTensorAlgebra`: an algebra on tuples of full Grassmannian RC graphs (as `CrystalGraphTensor`
keys), with conversion to `RCGraphRing` by squash-multiplying the factors together.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorPrintingTerm"></a>

## GrassTensorPrintingTerm Objects

```python
class GrassTensorPrintingTerm(PrintingTerm)
```

Display symbol for a `GrassTensorAlgebra` basis tuple.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebraElement"></a>

## GrassTensorAlgebraElement Objects

```python
class GrassTensorAlgebraElement(CrystalGraphRingElement)
```

Element of GrassTensorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra"></a>

## GrassTensorAlgebra Objects

```python
class GrassTensorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key: CrystalGraphTensor | tuple) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.grove_wc_ring"></a>

# schubmult.rings.combinatorial.grove\_wc\_ring

`GroveWCGraphRing` / `DualGroveWCGraphRing`: `WCGraphRing` quotients modeling grove polynomials
(the K-theoretic analogue of forest polynomials). WC graphs are snapped to canonical grove
representatives; products factor through `BoundedWCFactorAlgebra`.

<a id="schubmult.rings.combinatorial.grove_wc_ring.GroveWCGraphRingElement"></a>

## GroveWCGraphRingElement Objects

```python
class GroveWCGraphRingElement(WCGraphRingElement)
```

Element of `GroveWCGraphRing`.

<a id="schubmult.rings.combinatorial.grove_wc_ring.DualGroveWCGraphRingElement"></a>

## DualGroveWCGraphRingElement Objects

```python
class DualGroveWCGraphRingElement(GroveWCGraphRingElement)
```

Element of `DualGroveWCGraphRing`.

<a id="schubmult.rings.combinatorial.grove_wc_ring.GroveWCGraphRing"></a>

## GroveWCGraphRing Objects

```python
class GroveWCGraphRing(WCGraphRing)
```

`WCGraphRing` snapped to grove-class representatives; see the module docstring.

<a id="schubmult.rings.combinatorial.grove_wc_ring.DualGroveWCGraphRing"></a>

## DualGroveWCGraphRing Objects

```python
class DualGroveWCGraphRing(WCGraphRing)
```

Dual-product variant of `GroveWCGraphRing`.

<a id="schubmult.rings.combinatorial.hw_rc_ring"></a>

# schubmult.rings.combinatorial.hw\_rc\_ring

`HWRCGraphRing`: `RCGraphRing` quotient onto crystal highest-weight RC graphs.

Every product is snapped to its highest-weight representative, so basis elements
index Demazure crystal components. ``coproduct_on_basis`` is defined on highest-weight
elements via the free-algebra Schubert coproduct and `GrassTensorAlgebra`.

<a id="schubmult.rings.combinatorial.hw_rc_ring.HWRCGraphRing"></a>

## HWRCGraphRing Objects

```python
class HWRCGraphRing(RCGraphRing)
```

`RCGraphRing` with products projected onto highest-weight RC graphs; see the module docstring.

<a id="schubmult.rings.combinatorial.key_rc_ring"></a>

# schubmult.rings.combinatorial.key\_rc\_ring

`KeyRCGraphRing`: `RCGraphRing` quotient modeling key polynomials (Demazure characters).

RC graphs are snapped to a canonical representative of their ``extremal_weight`` class
(matching Edelman-Greene recording tableaux), and ``to_free_algebra_element`` lands in
the dual key basis indexed by ``extremal_weight``.

<a id="schubmult.rings.combinatorial.key_rc_ring.KeyRCGraphRingElement"></a>

## KeyRCGraphRingElement Objects

```python
class KeyRCGraphRingElement(RCGraphRingElement)
```

Element of `KeyRCGraphRing`; converts to the free-algebra key basis via ``extremal_weight``.

<a id="schubmult.rings.combinatorial.key_rc_ring.KeyRCGraphRingElement.to_free_algebra_element"></a>

#### to\_free\_algebra\_element

```python
def to_free_algebra_element(basis=None)
```

Map each RC graph to the dual key basis element indexed by its ``extremal_weight``.

<a id="schubmult.rings.combinatorial.key_rc_ring.KeyRCGraphRing"></a>

## KeyRCGraphRing Objects

```python
class KeyRCGraphRing(RCGraphRing)
```

`RCGraphRing` with products snapped to canonical key-class representatives; see the module docstring.

<a id="schubmult.rings.combinatorial.plactic_algebra"></a>

# schubmult.rings.combinatorial.plactic\_algebra

`PlacticAlgebra` and `NilPlacticAlgebra`: rings whose basis elements are `Plactic` /
`NilPlactic` tableaux (plactic and nilplactic monoid algebras).

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticPrintingTerm"></a>

## PlacticPrintingTerm Objects

```python
class PlacticPrintingTerm(TypedPrintingTerm)
```

Display symbol for a `PlacticAlgebra` basis tableau.

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticAlgebraElement"></a>

## PlacticAlgebraElement Objects

```python
class PlacticAlgebraElement(BaseRingElement)
```

PlacticAlgebra elements are linear combinations of Plactic basis elements.

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticAlgebra"></a>

## PlacticAlgebra Objects

```python
class PlacticAlgebra(BaseRing)
```

The plactic monoid algebra on `Plactic` tableaux (``op=True`` for the opposite product).

<a id="schubmult.rings.combinatorial.qy_rc_graph_ring"></a>

# schubmult.rings.combinatorial.qy\_rc\_graph\_ring

`QYRCGraphRing`: `RCGraphRing` quotient onto quasi-Yamanouchi RC graphs.

Every product is snapped by merging mergeable adjacent rows (``_canonical_rc``, the
same normalization as `RCGraph.snap_qy`), so basis elements are quasi-Yamanouchi.

<a id="schubmult.rings.combinatorial.qy_rc_graph_ring.QYRCGraphRing"></a>

## QYRCGraphRing Objects

```python
class QYRCGraphRing(RCGraphRing)
```

`RCGraphRing` with products snapped to quasi-Yamanouchi form; see the module docstring.

<a id="schubmult.rings.combinatorial.rc_graph_ring"></a>

# schubmult.rings.combinatorial.rc\_graph\_ring

`RCGraphRing`: the ring whose basis elements are `RCGraph`s.

Two products live here. ``*`` is the dual (stacking) product from `RCGraph.product`,
defined for every pair and compatible with the Schubert-polynomial coproduct
(``vertical_coproduct``). ``%`` (``rc_product``) is the polynomial product,
currently implemented only when one factor is Grassmannian with a sufficiently
large descent (or in two-row cases via squash decomposition). Crystal operators
and most `RCGraph` methods extend linearly to ring elements (see ``broadcast``).
``schub(perm, n)`` is the sum of all RC graphs of ``perm`` with ``n`` rows, i.e. the
Schubert polynomial as a ring element.

`GrassRCGraphRing` restricts to Grassmannian RC graphs (single descent at the last row).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement"></a>

## RCGraphRingElement Objects

```python
class RCGraphRingElement(CrystalGraphRingElement, SchubertMonomialRingElement)
```

RCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for %; the approach taken is to define the ambiguous term in the Leibniz formula
rather than compute the polynomial product directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.as_terms"></a>

#### as\_terms

```python
def as_terms()
```

Terms ``coeff * rc`` in dict order; the empty RC graph is kept as a symbol rather than collapsing to 1.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by RC graph (sympy printing hook).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.vex"></a>

#### vex

```python
@property
def vex()
```

Linear extension of `RCGraph.vex`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.grass"></a>

#### grass

```python
@property
def grass()
```

Linear extension of `RCGraph.grass`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.squash_product"></a>

#### squash\_product

```python
def squash_product(other_rc)
```

Linear extension of RCGraph.squash_product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.vertical_coproduct"></a>

#### vertical\_coproduct

```python
def vertical_coproduct()
```

Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
and induces the mul product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Coproduct via `RCGraphRing.coproduct_on_basis` (currently pattern-restricted; see there).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an RCGraphRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.almosteq"></a>

#### almosteq

```python
def almosteq(other)
```

Whether ``self - other`` has all-zero coefficients.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an RCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an RCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

Linear extension of RCGraph.crystal_reflection:
For each basis RCGraph, apply its crystal_reflection(index) and collect results.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.weight_reflection"></a>

#### weight\_reflection

```python
def weight_reflection(index)
```

Linear extension of `RCGraph.weight_reflection`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.shiftup"></a>

#### shiftup

```python
def shiftup(k)
```

Linear extension of `RCGraph.shiftup`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.prepend"></a>

#### prepend

```python
def prepend(k)
```

Linear extension of `RCGraph.prepend`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
def zero_out_last_row()
```

Linear extension of `RCGraph.zero_out_last_row`, dropping terms whose last row is nonempty.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.resize"></a>

#### resize

```python
def resize(n)
```

Linear extension of `RCGraph.resize`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.clip"></a>

#### clip

```python
def clip(n)
```

Keep the first ``n`` rows of each RC graph (left factor of `RCGraph.vertical_cut`).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.transpose"></a>

#### transpose

```python
def transpose(length)
```

Linear extension of `RCGraph.transpose`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.project"></a>

#### project

```python
def project()
```

Round-trip through the free algebra: ``from_free_algebra_element(to_free_algebra_element())``
(projects onto the sum-of-all-RC-graphs representatives).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.trim_operator"></a>

#### trim\_operator

```python
def trim_operator(i)
```

Linear extension of `RCGraphRing.trim_operator`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.double_elem_sym_squash"></a>

#### double\_elem\_sym\_squash

```python
def double_elem_sym_squash(weight, yvars, zvars)
```

Linear extension of `RCGraph.double_elem_sym_squash`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.full_double_elem_sym_squash"></a>

#### full\_double\_elem\_sym\_squash

```python
def full_double_elem_sym_squash(p, yvars, zvars)
```

Linear extension of `RCGraph.full_double_elem_sym_squash`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.grass_coaction"></a>

#### grass\_coaction

```python
def grass_coaction()
```

Linear extension of `RCGraphRing.grass_coaction`.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.broadcast"></a>

#### broadcast

```python
@property
def broadcast()
```

Proxy that lifts any `RCGraph` method linearly: ``elem.broadcast.method(*args)`` applies
``rc.method(*args)`` to each basis RC graph and sums the results with coefficients.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing"></a>

## RCGraphRing Objects

```python
class RCGraphRing(SchubertMonomialRing, CrystalGraphRing)
```

The ring of `RCGraph`s; see the module docstring. Instances are distinct (hashed by an id counter).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

Wrap an `RCGraph` as a basis element (or re-parent an existing element).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.new"></a>

#### new

```python
def new(x)
```

The basis element for RC graph ``x`` with coefficient 1.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.monomial"></a>

#### monomial

```python
def monomial(*tup)
```

The ring element for the monomial ``x^tup``: product of one-row RC graphs, split recursively.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.elem_sym"></a>

#### elem\_sym

```python
def elem_sym(descent, weight)
```

Sum of all RC graphs for the elementary-symmetric permutation with the given ``weight`` and descent.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.from_free_algebra_element"></a>

#### from\_free\_algebra\_element

```python
def from_free_algebra_element(elem)
```

Convert a free-algebra element (via the word basis) into a sum of ``monomial`` elements.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the RCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.double_schub"></a>

#### double\_schub

```python
def double_schub(perm, coeff_genset, n=None)
```

Return the RCGraphRing element corresponding to the double Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used), with `coeff_genset`
as the set of variables for the double part.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.trim_operator"></a>

#### trim\_operator

```python
def trim_operator(i, rc)
```

Divided difference at ``i`` followed by pulling out row ``i`` (only for terms where that row emptied).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.grass_coaction"></a>

#### grass\_coaction

```python
def grass_coaction(elem: RCGraph)
```

Coaction of the Grassmannian ring: squash-decompose ``elem`` as ``base * grass``, coproduct
``grass`` in `GrassRCGraphRing`, and reattach ``base`` to the left factors.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(rc)
```

Coproduct of a single RC graph, lifted from the free-algebra Schubert coproduct via
`BoundedRCFactorAlgebra`; only implemented for permutations avoiding ``4132``, ``1432``, ``3142``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.old_coproduct_on_basis"></a>

#### old\_coproduct\_on\_basis

```python
def old_coproduct_on_basis(elem)
```

Earlier recursive coproduct (peel the last row, recurse, correct); superseded by ``coproduct_on_basis``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.rc_product"></a>

#### rc\_product

```python
def rc_product(elem1, elem2)
```

Polynomial product (``%``), bilinear extension of ``rc_single_product``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.rc_single_product"></a>

#### rc\_single\_product

```python
def rc_single_product(u_rc, v_rc)
```

Polynomial product of two RC graphs of equal length: handled via squash decomposition for two
rows, and via ``squash_product``/``left_squash`` when one factor is Grassmannian with a large
enough descent; raises ``NotImplementedError`` otherwise.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.potential_products"></a>

#### potential\_products

```python
@cache
def potential_products(left, right, length)
```

Candidate RC graphs that could appear in ``left % right``, obtained by multiplying the
transposes with ``*`` and transposing back.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.potential_prodperms"></a>

#### potential\_prodperms

```python
def potential_prodperms(left, right, length)
```

Permutations of the ``potential_products``.

<a id="schubmult.rings.combinatorial.rc_graph_ring.GrassRCGraphRing"></a>

## GrassRCGraphRing Objects

```python
class GrassRCGraphRing(RCGraphRing)
```

`RCGraphRing` restricted to Grassmannian RC graphs (``_is_valid_grass``); non-Grassmannian terms
are filtered out of sums and products, and ``coproduct_on_basis`` is computed recursively by
vertical cuts.

<a id="schubmult.rings.combinatorial.rc_graph_ring.GrassRCGraphRing.mul_pair"></a>

#### mul\_pair

```python
def mul_pair(a, b)
```

Product of two Grassmannian RC graphs: stack ``b`` (shifted) below ``a`` and keep valid Grassmannian results.

<a id="schubmult.rings.combinatorial.rc_graph_ring.GrassRCGraphRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(elem)
```

Coproduct of a Grassmannian RC graph: one-row case is the standard Pieri splitting; otherwise
cut in half vertically, coproduct each half, and keep pairs whose ``%`` product recovers ``elem``.

<a id="schubmult.rings.combinatorial.rc_schubert_ring"></a>

# schubmult.rings.combinatorial.rc\_schubert\_ring

Orphaned earlier variant of `bounded_rc_factor_algebra` (defines the same class names, is not
imported anywhere, and is not exported from the package). Kept for reference; prefer
`schubmult.rings.combinatorial.bounded_rc_factor_algebra`.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorPrintingTerm"></a>

## BoundedRCFactorPrintingTerm Objects

```python
class BoundedRCFactorPrintingTerm(PrintingTerm)
```

Display symbol for a basis key (orphaned variant).

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebraElement"></a>

## BoundedRCFactorAlgebraElement Objects

```python
class BoundedRCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebraElement.prune"></a>

#### prune

```python
def prune()
```

Merge terms whose evaluated RC graphs share forest/omega invariants.

Terms are grouped by
``(rc.forest_weight, rc.omega_invariant[1])`` where
``rc = self.ring.key_to_rc_graph(key)``.
One key per group is kept (first encountered), and coefficients are summed.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra"></a>

## BoundedRCFactorAlgebra Objects

```python
class BoundedRCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring"></a>

# schubmult.rings.combinatorial.schubert\_monomial\_ring

Schubert Monomial Ring module

Provides base classes for rings whose basis elements represent Schubert monomials
(e.g., RC-graphs, BPDs, pipe dreams) with common operations like expansion to
polynomials, divided differences, and crystal operations.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialPrintingTerm"></a>

## SchubertMonomialPrintingTerm Objects

```python
class SchubertMonomialPrintingTerm(TypedPrintingTerm)
```

Printing term for Schubert monomial basis elements.

Delegates printing to the underlying key object (typically an RCGraph, BPD, etc.)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement"></a>

## SchubertMonomialRingElement Objects

```python
class SchubertMonomialRingElement(BaseRingElement)
```

Base class for ring elements whose basis elements are Schubert monomials.

This provides a common interface for objects like:
- RCGraphRingElement (basis elements are RCGraphs)
- BPDRingElement (basis elements are BPDs)

Common operations include:
- Polynomial expansion via polyvalue()
- Divided difference operators
- Crystal structure operations (if the basis elements support them)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x, y=None, *args, **kwargs)
```

Evaluate as a polynomial in variables x (and optionally y).

Linear extension: for each basis element, call its polyvalue() method
and sum the results weighted by coefficients.

**Arguments**:

- `x` - Variable or sequence of variables for polynomial evaluation
- `y` - Optional second set of variables for double Schubert polynomials
- ```**kwargs``` - Additional arguments passed to basis element polyvalue
  

**Returns**:

  Symbolic expression representing the polynomial

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms ``coeff * basis_symbol`` in dict order (sympy printing hook).

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.to_free_algebra_element"></a>

#### to\_free\_algebra\_element

```python
def to_free_algebra_element(basis=None, *, word=False)
```

Convert to FreeAlgebra element in Schubert basis.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing"></a>

## SchubertMonomialRing Objects

```python
class SchubertMonomialRing(BaseRing)
```

Base class for rings whose basis elements are Schubert monomials.

Inherits from BaseRing to provide standard ring operations (add, sub, mul, etc.)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.printing_term"></a>

#### printing\_term

```python
def printing_term(key)
```

Wrap the basis key in a `SchubertMonomialPrintingTerm`.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.from_dict"></a>

#### from\_dict

```python
def from_dict(dct)
```

Build an element from ``{key: coeff}`` without coefficient coercion.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.mul"></a>

#### mul

```python
def mul(a, b)
```

Multiply two elements via each basis key's ``product`` method (which returns ``{key: coeff}``),
or scale by a scalar ``b``.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.rmul"></a>

#### rmul

```python
def rmul(a, b)
```

Scale by the scalar ``b``.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(key)
```

The basis element for ``key`` with coefficient 1.

<a id="schubmult.rings.combinatorial.schubert_rc_ring"></a>

# schubmult.rings.combinatorial.schubert\_rc\_ring

`SchubertRCGraphRing`: `RCGraphRing` whose product is computed through `BoundedRCFactorAlgebra`.

Each RC graph is factored into elementary-symmetric RC graphs (``_factor_rc``), the factors
are multiplied in the bounded factor algebra, and the result is converted back; ``schubert_poly``
is the sum of all RC graphs of a permutation.

<a id="schubmult.rings.combinatorial.schubert_rc_ring.SchubertRCGraphRingElement"></a>

## SchubertRCGraphRingElement Objects

```python
class SchubertRCGraphRingElement(RCGraphRingElement)
```

Element of `SchubertRCGraphRing`.

<a id="schubmult.rings.combinatorial.schubert_rc_ring.SchubertRCGraphRing"></a>

## SchubertRCGraphRing Objects

```python
class SchubertRCGraphRing(RCGraphRing)
```

`RCGraphRing` multiplying via `BoundedRCFactorAlgebra` factorizations; see the module docstring.

<a id="schubmult.rings.combinatorial.slide_rc_ring"></a>

# schubmult.rings.combinatorial.slide\_rc\_ring

`SlideRCGraphRing`: `RCGraphRing` modeling fundamental slide polynomials.

RC graphs are snapped to a canonical representative of their quasi-Yamanouchi class
(``snap_qy().length_vector``), and ``slide_poly(comp)`` is the sum of RC graphs with
quasi-Yamanouchi weight ``comp``. Products go through `BoundedRCFactorAlgebra`.

<a id="schubmult.rings.combinatorial.slide_rc_ring.SlideRCGraphRingElement"></a>

## SlideRCGraphRingElement Objects

```python
class SlideRCGraphRingElement(RCGraphRingElement)
```

Element of `SlideRCGraphRing`.

<a id="schubmult.rings.combinatorial.slide_rc_ring.SlideRCGraphRing"></a>

## SlideRCGraphRing Objects

```python
class SlideRCGraphRing(RCGraphRing)
```

`RCGraphRing` snapped to quasi-Yamanouchi-class representatives; see the module docstring.

<a id="schubmult.rings.combinatorial.wc_graph_ring"></a>

# schubmult.rings.combinatorial.wc\_graph\_ring

`WCGraphRing`: the ring whose basis elements are `WCGraph`s (the K-theoretic / Grothendieck
analogue of `RCGraphRing`). ``to_free_algebra_element`` lands in the free-algebra Grothendieck
basis by default.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRingElement"></a>

## WCGraphRingElement Objects

```python
class WCGraphRingElement(SchubertMonomialRingElement)
```

WCGraphRing elements are linear combinations of WCGraph basis elements.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing"></a>

## WCGraphRing Objects

```python
class WCGraphRing(SchubertMonomialRing)
```

The ring of `WCGraph`s; see the module docstring.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing.groth"></a>

#### groth

```python
def groth(perm, n=None)
```

Return the WCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.direct_product_ring"></a>

# schubmult.rings.direct\_product\_ring

`DirectProductRing`: the direct product ``R_0 x R_1 x ... x R_n`` of `BaseRing` instances.

Keys are ``(i, k)`` with ``i`` the component index and ``k`` a key of ``R_i``; multiplication is
componentwise and cross-component products vanish. Use ``from_component``/``project`` to move
elements in and out, and ``elem[i]`` to read component ``i``.

<a id="schubmult.rings.direct_product_ring.DirectProductRing"></a>

## DirectProductRing Objects

```python
class DirectProductRing(BaseRing)
```

Direct product of an arbitrary (fixed) number of rings.

An element is stored as a dict mapping keys ``(i, k)`` to coefficients,
where ``i`` is the component index and ``k`` is a basis key from the
*i*-th ring.  Addition and multiplication are componentwise: terms from
different components never interact, and cross-component products are
zero.

Parameters
----------
``*rings`` : BaseRing
    One or more constituent rings.

Examples
--------
>>> D = DirectProductRing(R1, R2, R3)
>>> a = D.from_component(0, some_R1_element)
>>> b = D.from_component(2, some_R3_element)
>>> a + b          # lives in components 0 and 2
>>> a * b          # zero (different components)
>>> D[0]           # R1
>>> D.project(a, 0)  # back to R1

<a id="schubmult.rings.direct_product_ring.DirectProductRing.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(i)
```

Return the *i*-th constituent ring.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.rings"></a>

#### rings

```python
@property
def rings()
```

The tuple of component rings.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.one"></a>

#### one

```python
@property
def one()
```

The identity: the sum of every component's identity.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.component_one"></a>

#### component\_one

```python
def component_one(i)
```

Return the identity element supported only on component *i*.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.from_component"></a>

#### from\_component

```python
def from_component(i, elem)
```

Lift an element of ``self[i]`` into the direct product.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.project"></a>

#### project

```python
def project(elem, i)
```

Project onto component *i*, returning an element of ``self[i]``.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Componentwise multiplication.

Only terms in the same component interact; cross-component products
are zero.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(*elems)
```

Construct an element from one element per component.

``D(e0, e1, ..., en)`` lifts each ``ei`` (an element of ``D[i]``)
into the direct product and sums them.

<a id="schubmult.rings.direct_product_ring.DirectProductBasisElement"></a>

## DirectProductBasisElement Objects

```python
class DirectProductBasisElement(PrintingTerm)
```

Printing term for a key ``(i, k)``; renders as ``(term)_i``.

<a id="schubmult.rings.direct_product_ring.DirectProductRingElement"></a>

## DirectProductRingElement Objects

```python
class DirectProductRingElement(BaseRingElement)
```

Element of a `DirectProductRing`; ``elem[i]`` projects onto component ``i``.

<a id="schubmult.rings.direct_product_ring.DirectProductRingElement.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key)
```

Index by component integer or by ``(i, k)`` basis key.

* ``elem[i]`` — project onto component *i* (returns a ``self.ring[i]`` element).
* ``elem[(i, k)]`` — coefficient lookup (standard dict behavior).

<a id="schubmult.rings.free_algebra"></a>

# schubmult.rings.free\_algebra

The free algebra: the graded dual of `schubmult.rings.polynomial_algebra`.

The free (noncommutative) algebra on generators indexed by nonnegative integers is
the graded dual of the polynomial ring ``Z[x_1, x_2, ...]``, under the pairing in
which the word ``(a_1, ..., a_n)`` is dual to the monomial
``x_1^{a_1} x_2^{a_2} ... x_n^{a_n}`` (i.e. a word is the exponent vector of its dual
monomial). Concretely, ``FreeAlgebraElement.pairing(poly_elem)`` /
``PolynomialAlgebraElement.apply_dual_element(free_elem)`` sum the products of
matching coefficients once both sides are in the word / monomial basis.

Under this duality:

- concatenation of words (the `WordBasis` product) is dual to the variable-splitting
  coproduct on polynomials (`PolynomialAlgebraElement.branch`/``coproduct``), and the
  free algebra's coproduct is dual to polynomial multiplication;
- every free-algebra basis is the dual of a polynomial-algebra basis, exposed via
  ``Basis.dual_basis()``: `SchubertBasis` <-> ``SchubertPolyBasis``, `KeyBasis` <->
  ``KeyPolyBasis``, `ForestBasis` <-> ``ForestPolyBasis``, `FundamentalSlideBasis` <->
  ``FundamentalSlidePolyBasis``, `GrothendieckBasis` <-> ``GrothendieckPolyBasis``, and
  so on, with `WordBasis` <-> ``MonomialBasis`` as the pair everything is computed through.

Basis keys carry a *length* (number of variables) alongside the combinatorial index,
e.g. a `SchubertBasis` key is ``(perm, numvars)``: the dual of ``S_perm`` regarded as
a polynomial in exactly ``numvars`` variables. This grading by number of variables is
what makes the duality with words of a fixed length work.

The core classes are `FreeAlgebra` (the ring, parametrized by a `FreeAlgebraBasis`
class) and `FreeAlgebraElement`; ``change_basis`` converts between bases.

Pre-built instances:
    - ``FA``: `WordBasis` (the default)
    - ``ASx``: `SchubertBasis`
    - ``AGx``: `GrothendieckBasis`
    - ``ADSx``: the double Schubert separated-descents ring used for expansion
    - ``ForestDual``, ``GroveDual``, ``GlideDual``: the corresponding bases

Available bases:
    WordBasis, SchubertBasis, CompositionSchubertBasis, ElementaryBasis,
    ForestBasis, FundamentalSlideBasis, GlideBasis, GrothendieckBasis, GroveBasis,
    JBasis, JTBasis, KeyBasis, LascouxBasis, MonomialSlideBasis, NElementaryBasis,
    SchubertSchurBasis, SchurElementaryBasis, SeparatedDescentsBasis, ZBasis.

<a id="schubmult.rings.free_algebra._core"></a>

# schubmult.rings.free\_algebra.\_core

`FreeAlgebra` and `FreeAlgebraElement`: the graded dual of the polynomial algebra.

See the package docstring (`schubmult.rings.free_algebra`) for the duality. This module
holds the ring and element classes; the individual bases live in sibling modules and
plug in through the `FreeAlgebraBasis` interface. ``FA``, ``ASx``, ``AGx`` are the
standard instances.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement"></a>

## FreeAlgebraElement Objects

```python
class FreeAlgebraElement(BaseRingElement)
```

An element of a `FreeAlgebra`: ``{basis_key: coefficient}``.

In the `WordBasis` a key is a word -- a tuple of nonnegative integers -- dual to
the monomial with that exponent vector. In other bases the key is that basis's
combinatorial index together with a number of variables (e.g. ``(perm, numvars)``
for `SchubertBasis`). Beyond ring arithmetic, elements support basis changes
(``change_basis``), the duality pairing with polynomials (``pairing``,
``poly_inner_product``), expansion into Schubert rings (``schub_expand``), and
word-level operations (``inject``, ``prefix``, ``suffix``, ``interval``, ``split``,
``factorize``) that are computed in the word basis and transported back.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.interleave"></a>

#### interleave

```python
def interleave(other, zero_pad=True)
```

Interleave two elements letter-by-letter in the word basis.

Converts both elements to WordBasis, then interleaves each pair of
words by alternating entries (a1, b1, a2, b2, ...). Shorter words
are zero-padded when ``zero_pad`` is True.

**Arguments**:

- `other` - Another FreeAlgebraElement to interleave with.
- `zero_pad` - If True, pad shorter words with zeros.
  

**Returns**:

  The interleaved element in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.inject"></a>

#### inject

```python
def inject(i, other)
```

Insert another element's words at position *i* in this element's words.

Delegates to the current basis's ``inject`` classmethod.

**Arguments**:

- `i` - Nonnegative integer insertion index.
- `other` - Another FreeAlgebraElement to inject.
  

**Returns**:

  A new element with the injected words.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.prefix"></a>

#### prefix

```python
def prefix(length)
```

Extract the first *length* letters of each word.

Delegates to the current basis's ``prefix`` classmethod.

**Arguments**:

- `length` - Nonnegative integer prefix length.
  

**Returns**:

  A new element containing the prefixes.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.suffix"></a>

#### suffix

```python
def suffix(length)
```

Extract the last *length* letters of each word.

Delegates to the current basis's ``suffix`` classmethod.

**Arguments**:

- `length` - Nonnegative integer suffix length.
  

**Returns**:

  A new element containing the suffixes.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.interval"></a>

#### interval

```python
def interval(start, stop)
```

Extract a subword from position *start* to *stop* in each word.

Delegates to the current basis's ``interval`` classmethod.

**Arguments**:

- `start` - Nonnegative start index (inclusive).
- `stop` - Stop index (exclusive).
  

**Returns**:

  A new element containing the subwords.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.poly_inner_product"></a>

#### poly\_inner\_product

```python
def poly_inner_product(poly, genset, n)
```

The duality pairing of this element with a polynomial expression.

Expands ``poly`` into monomials in ``genset`` (exponent vectors padded/truncated to
``n`` variables, or trailing zeros stripped if ``n`` is ``None``), converts ``self``
to the word basis, and sums ``coeff_word * coeff_monomial`` over matching
word/exponent-vector pairs.

**Arguments**:

- `poly` - A polynomial expression.
- `genset` - The generating set of variables for the polynomial.
- `n` - Number of variables to use (or None for automatic).
  

**Returns**:

  The integer inner product value.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.kill_zero"></a>

#### kill\_zero

```python
def kill_zero(fat=False, val=S.Zero)
```

Remove zeros from each word key.

In the word basis, strips all zero entries from each key. When
*fat* is True, multiplies the coefficient by ``val`` raised to
the number of removed zeros instead of simply dropping them.

**Arguments**:

- `fat` - If True, weight by ``val`` per removed zero.
- `val` - The value to raise per zero when *fat* is True.
  

**Returns**:

  A new element with zeros removed, in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

Return a dict mapping printing terms to sympified coefficients.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand all coefficients symbolically.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.pairing"></a>

#### pairing

```python
def pairing(other)
```

The duality pairing with a `PolynomialAlgebraElement`.

Converts ``self`` to the word basis and ``other`` to the monomial basis, then sums
``coeff_word * coeff_monomial`` over words equal to exponent vectors. This is the
pairing under which the free algebra is the graded dual of the polynomial algebra.

**Arguments**:

- `other` - A polynomial algebra element to pair with.
  

**Returns**:

  The integer pairing value.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.hom_nsym"></a>

#### hom\_nsym

```python
def hom_nsym()
```

Apply the homomorphism to noncommutative symmetric functions.

Strips zero entries from each word key (dropping them from the word)
and returns the result in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_double_expand"></a>

#### tup\_double\_expand

```python
@staticmethod
@cache
def tup_double_expand(tup)
```

Realize the word ``tup`` as a double Schubert (separated-descents) ring element: the
product ``prod_i S_{uncode([tup[i]])}`` with each factor in its own variable.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Realize the word ``tup`` as a single Schubert (separated-descents) ring element: the
product ``prod_i S_{uncode([tup[i]])}`` with each factor in its own variable, computed
by divide-and-conquer. This is the map word -> ``h_{a_1}(x_1) h_{a_2}(x_2) ...`` that
sends the word basis onto complete-symmetric-in-one-variable products.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis)
```

Re-express this element in another `FreeAlgebraBasis`.

Uses ``self.ring._basis.transition(other_basis)``, which most bases implement by
routing through the `WordBasis`.

**Arguments**:

- `other_basis` - The target basis class (e.g. WordBasis, SchubertBasis).
  

**Returns**:

  A new FreeAlgebraElement in the target basis's ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_expand"></a>

#### schub\_expand

```python
def schub_expand()
```

Realize this element in the single Schubert separated-descents ring via ``tup_expand``
(each word becomes a product of one-variable complete symmetric functions).

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_double_expand"></a>

#### schub\_double\_expand

```python
def schub_double_expand()
```

Double-alphabet analogue of ``schub_expand`` via ``tup_double_expand``.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.bcoproduct"></a>

#### bcoproduct

```python
def bcoproduct()
```

The "bar" coproduct (see `WordBasis.bcoproduct`) in the tensor square ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.factorize"></a>

#### factorize

```python
def factorize(j)
```

Split each word at position *j*, returning a tensor element.

In the word basis, each word ``w`` maps to ``(w[:j], w[j:])``.
The result is expressed in the tensor ring of the original basis.

**Arguments**:

- `j` - Position at which to split each word.
  

**Returns**:

  An element of the tensor product ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.remove_zeros"></a>

#### remove\_zeros

```python
def remove_zeros(inserter=S.One)
```

Remove zero entries from each key, weighting by *inserter* per zero removed.

**Arguments**:

- `inserter` - Scalar multiplied per removed zero (default 1).
  

**Returns**:

  A new element with zeros stripped from keys.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.__truediv__"></a>

#### \_\_truediv\_\_

```python
def __truediv__(other)
```

Skew by a permutation: ``elem / u`` applies ``skew_element(w, u, n)`` to each ``(w, n)`` key
(the dual of multiplying by ``S_u`` on the polynomial side).

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.split"></a>

#### split

```python
def split(p)
```

Split each word at position *p* into a tensor element.

Words shorter than *p* are placed entirely in the left factor.

**Arguments**:

- `p` - Position at which to split.
  

**Returns**:

  An element of the tensor product ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.to_schub"></a>

#### to\_schub

```python
def to_schub(sym=False)
```

Convert to a Schubert ring element via ``tup_to_schub``.

**Arguments**:

- `sym` - If True, use symmetric expansion.
  

**Returns**:

  An element of the single Schubert ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra"></a>

## FreeAlgebra Objects

```python
class FreeAlgebra(BaseRing)
```

The free algebra on generators indexed by nonnegative integers, in a chosen basis.

The ring is basis-agnostic; a `FreeAlgebraBasis` *class* (not instance) supplies the
key type, product, coproduct, and transitions. ``FreeAlgebra(WordBasis)`` is the
concatenation algebra on words; ``FreeAlgebra(SchubertBasis)`` is the same algebra
written in the basis dual to Schubert polynomials. See the package docstring for
the duality with the polynomial algebra.

**Arguments**:

- `basis` - The basis class to use (default ``WordBasis``).
- `domain` - Coefficient domain (default ``EXRAW``).

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply every coefficient of *elem* by the scalar *x*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.tensor_schub_expand"></a>

#### tensor\_schub\_expand

```python
def tensor_schub_expand(tensor)
```

Expand a tensor element into the Schubert polynomial tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.__init__"></a>

#### \_\_init\_\_

```python
def __init__(basis=WordBasis, domain=None)
```

Initialize a FreeAlgebra with the given basis and coefficient domain.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.right_pad"></a>

#### right\_pad

```python
@staticmethod
def right_pad(tup, n)
```

Right-pad *tup* with zeros to length *n*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(key)
```

Compute the coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.bcoproduct_on_basis"></a>

#### bcoproduct\_on\_basis

```python
@cache
def bcoproduct_on_basis(key)
```

Compute the bar-coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via the basis product rule.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.from_rc_graph"></a>

#### from\_rc\_graph

```python
def from_rc_graph(rc_graph)
```

Create an element from an RC graph.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.matmul"></a>

#### matmul

```python
def matmul(elem, other)
```

Internal (Kronecker) product (``@`` operator) or scalar multiplication.

If ``other`` is a scalar, multiplies all coefficients. If it is a `FreeAlgebraElement`,
computes the basis's ``internal_product`` -- essentially the Kronecker product of
noncommutative symmetric functions, enumerated in the word basis by nonnegative
integer matrices with prescribed row/column sums (see `WordBasis.internal_product`;
requires SageMath).

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.new"></a>

#### new

```python
def new(*x)
```

Create a new element from the given key or arguments.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Return the display symbol for basis key *k*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.from_dict"></a>

#### from\_dict

```python
def from_dict(element)
```

Construct an element from a dict of ``{key: coefficient}`` pairs.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.skew_element"></a>

#### skew\_element

```python
def skew_element(w, u, n)
```

The skew element ``S_w / S_u`` in ``n`` variables (dual to multiplication by ``S_u``); see `SchubertBasis.skew_element`.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce a raw value into the coefficient domain.

<a id="schubmult.rings.free_algebra._core.make_AGx"></a>

#### make\_AGx

```python
def make_AGx(beta=None)
```

Create a Grothendieck-basis free algebra, optionally with custom beta.

<a id="schubmult.rings.free_algebra.composition_schubert_basis"></a>

# schubmult.rings.free\_algebra.composition\_schubert\_basis

`CompositionSchubertBasis`: `SchubertBasis` re-indexed by compositions. The key ``c`` (a code
padded with zeros to length ``numvars``) stands for the Schubert key ``(uncode(c), len(c))``,
so the dual polynomial basis is again Schubert polynomials.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis"></a>

## CompositionSchubertBasis Objects

```python
class CompositionSchubertBasis(FreeAlgebraBasis)
```

Schubert basis indexed by padded trimcode compositions.

A key is a composition ``c`` and corresponds to Schubert key
``(uncode(c), len(c))``.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list (composition).

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.as_schubert_key"></a>

#### as\_schubert\_key

```python
@classmethod
def as_schubert_key(cls, key)
```

Convert a composition key to a Schubert key ``(Permutation, length)``.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, key)
```

Normalize a key to a composition tuple.

Accepts either a Schubert key ``(Permutation, int)`` or a raw tuple.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the composition key for the given RC graph.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Inject *key2* into *key1* at position *i* using Schubert multiplication.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two composition keys by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
def coproduct(cls, key)
```

Compute the coproduct by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
def bcoproduct(cls, key)
```

Compute the bar-coproduct by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

Compute the skew element by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual basis (delegates to SchubertBasis).

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from CompositionSchubertBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``CompSchub``-labelled display object for the composition key *k*.

<a id="schubmult.rings.free_algebra.elementary_basis"></a>

# schubmult.rings.free\_algebra.elementary\_basis

`ElementaryBasis`: free-algebra basis indexed by ``(composition, numvars)``, dual to products of
elementary symmetric polynomials ``e_{c_1}(x_1..x_k) e_{c_2}(x_1..x_{k-1}) ...`` in nested
variable sets. `SchubertBasis` expands into it via the monomials of ``S_{perm * w0}``.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis"></a>

## ElementaryBasis Objects

```python
class ElementaryBasis(FreeAlgebraBasis)
```

Elementary symmetric function basis of the free algebra.

Keys are ``(tuple, int)`` pairs where the tuple encodes an elementary
symmetric function composition and the integer is the number of variables.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(tuple/list, int)`` pair.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, int)`` key.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, tup, numvars)
```

Transition an elementary key to the Schubert basis.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``Elem``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.forest_basis"></a>

# schubmult.rings.free\_algebra.forest\_basis

`ForestBasis`: the free-algebra basis dual to forest polynomials (``ForestPolyBasis``).

Keys are weak compositions (indexed-forest weights). Expansion into `SchubertBasis`
enumerates RC graphs and keeps those whose forest weight is the key. ``ForestDual`` is the
standard instance.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis"></a>

## ForestBasis Objects

```python
class ForestBasis(FreeAlgebraBasis)
```

Forest basis of the free algebra.

Keys are tuples representing indexed-forest weight vectors.
Transitions to the Schubert basis use RC graph enumeration
filtered by forest weight.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``ForestDual``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the ForestPolyBasis as the dual of ForestBasis.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a forest key to the Schubert basis via RC graph enumeration.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ForestBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.free_algebra_basis"></a>

# schubmult.rings.free\_algebra.free\_algebra\_basis

`FreeAlgebraBasis`: the interface a basis must implement to plug into `FreeAlgebra`.

A basis is a *class* (its methods are classmethods) defining a key type (``is_key``/
``as_key``/``zero_monom``), how to print a key, ``transition(other_basis)`` returning a
key -> ``{key: coeff}`` function into another basis, and ``dual_basis()`` naming the
`schubmult.rings.polynomial_algebra` basis it is dual to. Products, coproducts, and
the word-level operations all have default implementations that route through the
`WordBasis` via ``compose_transition``.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis"></a>

## FreeAlgebraBasis Objects

```python
class FreeAlgebraBasis()
```

Abstract base for free-algebra bases; see the module docstring.

Subclasses override the key methods (``is_key``, ``as_key``, ``zero_monom``,
``printing_term``, ``transition``, ``dual_basis``) and may override ``product``/
``coproduct`` with a direct rule; otherwise everything is computed in the `WordBasis`
and transported back.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid key for this basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Convert an RC graph to a basis-keyed dict.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a canonical key for this basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a function mapping keys of this basis to dicts in *other_basis*.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, key)
```

Return the display symbol for *key*.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.compose_transition"></a>

#### compose\_transition

```python
@classmethod
def compose_transition(cls, tkeyfunc, output)
```

Apply a key-level transition function to each key in *output*.

For each ``(key, v)`` in *output*, expands ``tkeyfunc(key)`` and
accumulates the results weighted by *v*.

**Arguments**:

- `tkeyfunc` - A function mapping a key to a ``{key: coeff}`` dict.
- `output` - A ``{key: coeff}`` dict to transform.
  

**Returns**:

  A merged ``{key: coeff}`` dict in the target basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

The `schubmult.rings.polynomial_algebra` basis this basis is dual to under the word/monomial pairing.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.change_tensor_basis"></a>

#### change\_tensor\_basis

```python
@classmethod
def change_tensor_basis(cls, tensor_elem, basis1, basis2)
```

Change the bases of both factors of a tensor element.

**Arguments**:

- `tensor_elem` - An element of a tensor product ring.
- `basis1` - Target basis for the left factor.
- `basis2` - Target basis for the right factor.
  

**Returns**:

  The tensor element re-expressed in the new bases.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of *key* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key)
```

Compute the bar-coproduct of *key* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

The internal (Kronecker) product of NSym (see `WordBasis.internal_product`), computed via the word basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Inject *key2* into *key1* at position *i* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.prefix"></a>

#### prefix

```python
@classmethod
def prefix(cls, key, length, coeff=S.One)
```

Extract a prefix of *length* letters by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.suffix"></a>

#### suffix

```python
@classmethod
def suffix(cls, key, length, coeff=S.One)
```

Extract a suffix of *length* letters by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.interval"></a>

#### interval

```python
@classmethod
def interval(cls, key, start, stop, coeff=S.One)
```

Extract a subword from *start* to *stop* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name)
```

Lazily resolve basis classes to avoid circular imports between basis modules.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis"></a>

# schubmult.rings.free\_algebra.fundamental\_slide\_basis

`FundamentalSlideBasis`: the free-algebra basis dual to fundamental slide polynomials
(``FundamentalSlidePolyBasis``). Keys are weak compositions; transitions are transposed from
the polynomial-side slide expansions.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis"></a>

## FundamentalSlideBasis Objects

```python
class FundamentalSlideBasis(FreeAlgebraBasis)
```

Fundamental slide basis of the free algebra.

Keys are weak composition tuples. Transitions are computed via
polynomial algebra slide polynomial expansions.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``FS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a fundamental slide key to the Schubert basis.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, key)
```

Transition a fundamental slide key to the word basis.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the FundamentalSlidePolyBasis as the dual.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from FundamentalSlideBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.glide_basis"></a>

# schubmult.rings.free\_algebra.glide\_basis

`GlideBasis`: the free-algebra basis dual to glide polynomials (``GlidePolyBasis``), the
K-theoretic analogue of `FundamentalSlideBasis`. Keys are weak compositions. ``GlideDual`` is
the standard instance.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis"></a>

## GlideBasis Objects

```python
class GlideBasis(FreeAlgebraBasis)
```

Glide basis of the free algebra.

Keys are weak composition tuples. Transitions are computed via
polynomial algebra slide polynomial expansions.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``FS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a glide key to the Grothendieck basis.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the GlidePolyBasis as the dual.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from GlideBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.grothendieck_basis"></a>

# schubmult.rings.free\_algebra.grothendieck\_basis

`GrothendieckBasis`: the free-algebra basis dual to Grothendieck polynomials.

Keys are ``(perm, numvars)`` as in `SchubertBasis`; the key is dual to ``G_perm`` in ``numvars``
variables (``GrothendieckPolyBasis`` on the polynomial side). Grothendieck polynomials are
taken at ``beta = 1``, which is without loss of generality: ``beta`` is recovered by grading
(a term of ``G_w`` in degree ``inv(w) + d`` carries ``beta^d``). The change of
basis to `SchubertBasis` is the transpose of the Grothendieck-to-Schubert expansion and is
computed combinatorially: enumerate unreduced BPDs of ``perm * w0``, take the co-BPD, and keep
the reduced ones, with sign ``(-1)^(inv(perm) - inv(result))``. Products and all other
transitions route through `SchubertBasis`. ``AGx`` is the standard instance.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis"></a>

## GrothendieckBasis Objects

```python
class GrothendieckBasis(FreeAlgebraBasis)
```

Free-algebra basis dual to Grothendieck polynomials (``beta = 1``); keys are ``(Permutation, numvars)``.

See the module docstring. Products and transitions go through `SchubertBasis`.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Whether ``x`` is ``(perm,)`` or ``(perm, numvars)``.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize to ``(Permutation, numvars)``; ``numvars`` defaults to the last descent.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
@cache
def transition_schubert(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in `SchubertBasis`.

For each unreduced BPD of ``perm * w0`` whose co-BPD is reduced with permutation ``u``
fitting in ``numvars`` variables, contributes ``(-1)^(inv(perm) - inv(u))`` to ``(u, numvars)``.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Key -> ``{key: coeff}`` function into ``other_basis``: identity on Grothendieck subclasses,
`transition_schubert` for `SchubertBasis`, and Schubert-then-onward for everything else.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Display as ``AGx(perm, numvars)``.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply by expanding both keys in `SchubertBasis`, using its separated-descents product,
and converting the result back.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

``GrothendieckPolyBasis``: Grothendieck polynomials are the dual basis.

<a id="schubmult.rings.free_algebra.grove_basis"></a>

# schubmult.rings.free\_algebra.grove\_basis

`GroveBasis`: the free-algebra basis dual to grove polynomials (``GrovePolyBasis``), the
K-theoretic analogue of `ForestBasis`. Keys are weak compositions (grove weights); expansion
into `GrothendieckBasis` enumerates WC graphs by grove weight. ``GroveDual`` is the standard
instance.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis"></a>

## GroveBasis Objects

```python
class GroveBasis(FreeAlgebraBasis)
```

Grove basis of the free algebra.

Keys are tuples representing indexed-grove weight vectors.
Transitions to the Grothendieck basis use RC graph enumeration
filtered by grove weight.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``GroveDual``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the GrovePolyBasis as the dual of GroveBasis.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a grove key to the Grothendieck basis via WC graph enumeration.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from GroveBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.j_basis"></a>

# schubmult.rings.free\_algebra.j\_basis

`JBasis`: free-algebra basis indexed by compositions with no zeros.

A Schubert key ``(perm, n)`` whose padded code has no zeros is itself a J key; zeros are
handled by the transitions in `SchubertBasis.transition_jbasis` and `WordBasis.transition_jbasis`.

<a id="schubmult.rings.free_algebra.j_basis.JBasis"></a>

## JBasis Objects

```python
class JBasis(FreeAlgebraBasis)
```

J basis of the free algebra.

Keys are tuples of positive integers (no zeros allowed in transitions).
The J basis indexes elements whose Schubert expansion has no zero
entries in the Lehmer code.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a J basis key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.coproduct"></a>

#### coproduct

```python
@classmethod
def coproduct(cls, key)
```

Coproduct for JBasis equals the bar-coproduct.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``J``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from JBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.jt_basis"></a>

# schubmult.rings.free\_algebra.jt\_basis

`JTBasis`: `JBasis` with a formal parameter ``t`` recording the number of stripped zeros.
Keys are ``(composition, power_of_t)``.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis"></a>

## JTBasis Objects

```python
class JTBasis(FreeAlgebraBasis)
```

JT basis of the free algebra (J basis with a parameter *t*).

Keys are ``(tuple, int)`` pairs where the tuple is a nonzero code
and the integer tracks a power of the parameter *t*.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a JT key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.pare_schubert"></a>

#### pare\_schubert

```python
@staticmethod
def pare_schubert(perm)
```

Extract a nonzero trimcode from *perm*, or None if it contains zeros.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.normalize_dct"></a>

#### normalize\_dct

```python
@staticmethod
def normalize_dct(dct)
```

Normalize a word dict by collecting zeros into leading positions.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a *t*-weighted ``JT``-labelled display object.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from JTBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.key_basis"></a>

# schubmult.rings.free\_algebra.key\_basis

`KeyBasis`: the free-algebra basis dual to key polynomials (Demazure characters).

Keys are weak compositions; the element is dual to the key polynomial with that weight
(``KeyPolyBasis``). Expansion into `SchubertBasis` enumerates RC graphs and keeps those whose
length vector is the extremal weight.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis"></a>

## KeyBasis Objects

```python
class KeyBasis(FreeAlgebraBasis)
```

Key polynomial (Demazure character) basis of the free algebra.

Keys are weak composition tuples. Transitions to the Schubert basis
use RC graph enumeration filtered by extremal weight.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the KeyPolyBasis as the dual of KeyBasis.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Key``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a key composition to the Schubert basis via RC graphs.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from KeyBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.lascoux_basis"></a>

# schubmult.rings.free\_algebra.lascoux\_basis

`LascouxBasis`: the free-algebra basis dual to Lascoux polynomials (``LascouxPolyBasis``),
the K-theoretic analogue of `KeyBasis`. Keys are weak compositions.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis"></a>

## LascouxBasis Objects

```python
class LascouxBasis(FreeAlgebraBasis)
```

Lascoux polynomial (Demazure character) basis of the free algebra.

Lascouxs are weak composition tuples. Transitions to the Schubert basis
use RC graph enumeration filtered by extremal weight.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the LascouxPolyBasis as the dual of LascouxBasis.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Lascoux``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a Lascoux composition to the Grothendieck basis via WC graphs.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from LascouxBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis"></a>

# schubmult.rings.free\_algebra.monomial\_slide\_basis

`MonomialSlideBasis`: the free-algebra basis dual to monomial slide polynomials. Keys are
weak compositions; transitions use coarsenings of compositions.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis"></a>

## MonomialSlideBasis Objects

```python
class MonomialSlideBasis(FreeAlgebraBasis)
```

Monomial slide basis of the free algebra.

Keys are weak composition tuples. Transitions use monomial slide
polynomial expansions and coarsenings of compositions.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``MS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
@classmethod
@cache
def transition_fundamental_slide(cls, key)
```

Transition a monomial slide key to the fundamental slide basis.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from MonomialSlideBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, key)
```

Transition a monomial slide key to the word basis.

<a id="schubmult.rings.free_algebra.nelementary_basis"></a>

# schubmult.rings.free\_algebra.nelementary\_basis

`NElementaryBasis`: the noncommutative elementary basis ``L`` of NSym inside the free algebra.

Keys are compositions (positive integers). ``L_alpha`` expands in words as the signed sum
``sum_{beta refines alpha} (-1)^(|alpha| - len(beta)) beta`` (refinements via SageMath), and the
product is concatenation.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis"></a>

## NElementaryBasis Objects

```python
class NElementaryBasis(FreeAlgebraBasis)
```

Non-commutative elementary basis (L basis) of the free algebra.

Keys are tuples of positive integers. The transition to the word
basis uses composition refinements from SageMath.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Concatenate two keys.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``L``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, tup)
```

Transition an NElementary key to the word basis via composition refinements (requires SageMath).

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from NElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schubert_basis"></a>

# schubmult.rings.free\_algebra.schubert\_basis

`SchubertBasis`: the free-algebra basis dual to Schubert polynomials.

A key is ``(perm, numvars)``: the element dual to ``S_perm`` viewed as a polynomial in
exactly ``numvars`` variables (so ``numvars >= max_descent(perm)``). Under the
word/monomial pairing this is the ``SchubertPolyBasis`` of the polynomial algebra.

The product is the separated-descents product (`SeparatedDescentsRing`): ``(u, p) * (v, q)``
places ``u`` in the first ``p`` variables and ``v`` in the next ``q``, giving a ``(w, p + q)``
expansion. ``transition_word`` expands a key into words via the SEM (elementary symmetric)
factorization of ``S_perm``, and the other ``transition_*`` methods reach the remaining
bases either directly or by way of the word basis. ``ASx`` is the standard instance.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis"></a>

## SchubertBasis Objects

```python
class SchubertBasis(FreeAlgebraBasis)
```

Free-algebra basis dual to Schubert polynomials; keys are ``(Permutation, numvars)``.

See the module docstring. Products go through the separated-descents Schubert ring,
and transitions to the word basis go through elementary symmetric function
decompositions.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Whether ``x`` is ``(perm,)`` or ``(perm, numvars)`` with ``perm`` a permutation/list/tuple.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

The key ``(rc_graph.perm, len(rc_graph))``: an RC graph's permutation with its row count as ``numvars``.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize to ``(Permutation, numvars)``; if ``numvars`` is omitted it defaults to the last descent.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Separated-descents product: ``(u, p) * (v, q)`` with ``u`` in the first ``p`` variables and
``v`` in the next ``q``, computed in `SeparatedDescentsRing`.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

The skew element ``S_w / S_u`` in ``n`` variables: the dual of multiplying by ``S_u``,
computed with the descent-side kernel ``schubmult_py_down`` and truncated to permutations
fitting in ``n`` variables.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Coproduct of ``(perm, numvars)`` (dual to polynomial multiplication): expand to words,
apply the word coproduct, and convert each tensor factor back to Schubert keys.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in the `GrothendieckBasis`, by taking the co-BPD of every RC graph
of ``perm * w0`` and collecting the resulting permutations.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schubert_schur"></a>

#### transition\_schubert\_schur

```python
@classmethod
def transition_schubert_schur(cls, *x)
```

Expand ``(perm, numvars)`` in the `SchubertSchurBasis`: split off the variables beyond
``numvars`` via a Schubert coproduct against a dominant permutation, yielding
``(partition, perm', numvars)`` keys.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schur_elementary"></a>

#### transition\_schur\_elementary

```python
@classmethod
def transition_schur_elementary(cls, *x)
```

Expand ``(perm, numvars)`` in the `SchurElementaryBasis` (a word-like tuple paired with a partition).

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_elementary"></a>

#### transition\_elementary

```python
@classmethod
def transition_elementary(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in the `ElementaryBasis`: read the monomials of ``S_{perm * w0}``
and complement each exponent against the staircase to get elementary-symmetric indices.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_separated_descents"></a>

#### transition\_separated\_descents

```python
@classmethod
def transition_separated_descents(cls, k, *x)
```

Expand ``(perm, numvars)`` in the level-``k`` `SeparatedDescentsBasis` via a Schubert coproduct
splitting the last ``k - 1`` variables, yielding ``(perm_left, perm_right, numvars)`` keys.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, perm, n)
```

Expand ``(perm, n)`` in the `JBasis`: a code with no zeros is already a J key; leading zeros are
peeled off (each contributing a factor ``t``, currently ``1``), and anything else goes via words.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

``SchubertPolyBasis``: Schubert polynomials are the dual basis under the word/monomial pairing.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition"></a>

#### transition

```python
@classmethod
@cache
def transition(cls, other_basis)
```

Return the key -> ``{key: coeff}`` function into ``other_basis``.

Direct routes exist for the word, elementary, Schubert-Schur, Schur-elementary,
composition-Schubert, separated-descents, and Grothendieck bases; everything else
is reached by going through the word basis first.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.old_transition_word"></a>

#### old\_transition\_word

```python
@classmethod
@cache
def old_transition_word(cls, perm, numvars)
```

Transition to WordBasis via SEM basis (legacy implementation).

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in the word basis.

Multiplies ``perm`` by the inverse of the dominant permutation with code
``(inv + numvars, ..., inv + 1)`` (``inv = perm.inv``) and writes the result in the SEM
basis with a custom ``elem_func`` that records each factor ``e_p(x_1..x_k)`` as the word
with ``k - p`` in position ``numvars - k + inv``. Products of such words in the monomial
polynomial algebra add letterwise, so the resulting polynomial *is* the word expansion.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Display symbol for ``(perm, numvars)`` (the separated-descents ``SepDescSchubPoly`` form).

<a id="schubmult.rings.free_algebra.schubert_schur_basis"></a>

# schubmult.rings.free\_algebra.schubert\_schur\_basis

`SchubertSchurBasis`: free-algebra basis dual to products ``s_lambda(x_1..x_n) * S_perm``
of a Schur polynomial in the first ``n`` variables with a Schubert polynomial. Keys are
``(partition, perm, numvars)``.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis"></a>

## SchubertSchurBasis Objects

```python
class SchubertSchurBasis(FreeAlgebraBasis)
```

Schubert-Schur basis of the free algebra.

Keys are ``(partition_tuple, Permutation)`` pairs encoding products of
a Grassmannian Schubert polynomial with a Schur polynomial.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(list/tuple, Permutation/list/tuple)`` pair.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, Permutation)`` key.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert-Schur key via the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, lambd, perm)
```

Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, lambd, perm)
```

Transition a Schubert-Schur key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from SchubertSchurBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``SS``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.schur_elementary_basis"></a>

# schubmult.rings.free\_algebra.schur\_elementary\_basis

`SchurElementaryBasis`: free-algebra basis dual to products of a nested elementary monomial
(as in `ElementaryBasis`) with a Schur polynomial. Keys are ``(elementary_tuple, partition)``.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis"></a>

## SchurElementaryBasis Objects

```python
class SchurElementaryBasis(FreeAlgebraBasis)
```

Schur-Elementary basis of the free algebra.

Keys are ``(tuple, tuple)`` pairs encoding products of
a standard elementary monomial (first tuple) and a Schur polynomial (second tuple, partition
of length precisely len(first_tuple) + 1 in increasing order, with zeros at the beginning if
needed).  For example, the key ``((1, 2), (0, 1, 3))`` corresponds to the product of the elementary monomial

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(list/tuple, list/tuple)`` pair.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, tuple)`` key.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert-Schur key via the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, elem_tup, lambd)
```

Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, elem_tup, lambd)
```

Transition a Schubert-Schur key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from SchurElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``SE``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.separated_descents_basis"></a>

# schubmult.rings.free\_algebra.separated\_descents\_basis

`SeparatedDescentsBasis(k)`: level-``k`` refinements of `SchubertBasis` in which a key
``(u, v, numvars)`` splits the Schubert index into a factor ``u`` and a factor ``v`` whose
descents are separated at ``k`` (the basis dual to the separated-descents factorization of
`schubmult.rings.schubert.separated_descents`).

`SeparatedDescentsBasis` is a factory producing a ``_SeparatedDescentsBasis`` subclass with
class attribute ``k``; products go through `SchubertBasis`, and `SchubertBasis` expands into
this basis via `SchubertBasis.transition_separated_descents`.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis"></a>

## \_SeparatedDescentsBasis Objects

```python
class _SeparatedDescentsBasis(FreeAlgebraBasis)
```

Separated descents basis of the free algebra (parameterized by level *k*).

Keys are ``(Permutation, Permutation, int)`` triples representing a
factorization into descents above and below a cutoff level *k*.
Instances are created by the :func:`SeparatedDescentsBasis` factory.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid separated descents key.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(Permutation, Permutation, int)`` key.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two separated descents keys via the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, perm0, perm1, numvars)
```

Transition a separated descents key to the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, perm0, perm1, n)
```

Transition a separated descents key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from this separated descents basis to *other_basis*.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``SepDesc<k>``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.separated_descents_basis.SeparatedDescentsBasis"></a>

#### SeparatedDescentsBasis

```python
def SeparatedDescentsBasis(k)
```

Factory that creates a separated descents basis class for level *k*.

<a id="schubmult.rings.free_algebra.word_basis"></a>

# schubmult.rings.free\_algebra.word\_basis

`WordBasis`: the word (concatenation) basis of the free algebra, dual to the monomial basis.

A key is a word ``(a_1, ..., a_n)`` of nonnegative integers, dual to the monomial
``x_1^{a_1} ... x_n^{a_n}``. The product is concatenation; the coproduct splits each
letter ``a`` into ``(i, a - i)`` (dual to polynomial multiplication). This is the hub
basis: every other `FreeAlgebraBasis` implements its operations by transitioning to
words and back, and this module holds the ``transition_*`` routines from words into
each of the other bases.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis"></a>

## WordBasis Objects

```python
class WordBasis(FreeAlgebraBasis)
```

Word basis of the free algebra: keys are tuples of nonnegative integers (words), each
dual to the monomial whose exponent vector is that word. Product is concatenation.
See the module docstring.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the length vector of the RC graph as a word key.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Concatenate two words (dual to the variable-splitting coproduct on polynomials).

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Insert *key2* into *key1* at position *i*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.prefix"></a>

#### prefix

```python
@classmethod
def prefix(cls, key, length, coeff=S.One)
```

Return the first *length* letters of *key*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.suffix"></a>

#### suffix

```python
@classmethod
def suffix(cls, key, length, coeff=S.One)
```

Return the last *length* letters of *key*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.interval"></a>

#### interval

```python
@classmethod
def interval(cls, key, start, stop, coeff=S.One)
```

Return the subword ``key[start:stop]``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key, coeff=S.One)
```

The coproduct of a word, dual to polynomial multiplication.

Each letter ``a`` splits into all ``(i, a - i)`` pairs (``x_j^a`` is the sum over
ways to write it as ``x_j^i * x_j^{a-i}``); the pieces are combined letterwise
by a divide-and-conquer tensor product.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key, coeff=S.One)
```

The "bar" coproduct: like `coproduct` but a zero letter is dropped rather than kept as
a ``0`` in the word, so word lengths are not preserved.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.try_internal_product"></a>

#### try\_internal\_product

```python
@classmethod
def try_internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product via integer matrices (requires SageMath).

Uses shifted keys (incremented by 1) with ``IntegerMatrices``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

The internal (Kronecker) product of two compositions (words without zeros), as in
noncommutative symmetric functions: sum over nonnegative integer matrices with row
sums ``key1`` and column sums ``key2`` of the word read off the nonzero entries.
Requires SageMath's ``IntegerMatrices``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a bracket-notation symbol like ``[210]`` for the word *k*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Expand a word tuple into the single Schubert basis via divide-and-conquer.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.jbasis_tup_expand"></a>

#### jbasis\_tup\_expand

```python
@staticmethod
@cache
def jbasis_tup_expand(tup)
```

Expand a word tuple into the Z basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a word key to the Schubert basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, key)
```

Transition a word key to the J basis via Pieri products.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_jtbasis"></a>

#### transition\_jtbasis

```python
@classmethod
def transition_jtbasis(cls, key)
```

Transition a word key to the JT basis via normalization.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_forest"></a>

#### transition\_forest

```python
@classmethod
def transition_forest(cls, key)
```

Transition a word key to the forest basis via RC graph enumeration.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_grove"></a>

#### transition\_grove

```python
@classmethod
def transition_grove(cls, key)
```

Transition a word key to the grove basis via WC graph enumeration.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the MonomialBasis as the dual of WordBasis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_monomial_slide"></a>

#### transition\_monomial\_slide

```python
@classmethod
def transition_monomial_slide(cls, key)
```

Transition a word key to the monomial slide basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_zbasis"></a>

#### transition\_zbasis

```python
@classmethod
def transition_zbasis(cls, key)
```

Expand a word in `ZBasis` by triangular elimination: repeatedly peel off the smallest
remaining word, subtracting the word expansion of the corresponding Z element.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_nelementary"></a>

#### transition\_nelementary

```python
@classmethod
def transition_nelementary(cls, tup)
```

Expand a composition in `NElementaryBasis`: signed sum over its refinements
(``(-1)^(|tup| - len(beta))``), via SageMath's ``Composition.finer``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_key"></a>

#### transition\_key

```python
@classmethod
def transition_key(cls, key)
```

Expand a word in `KeyBasis`: count RC graphs of length vector ``key`` whose extremal
weight equals their permutation's padded code (the dual of the key-to-monomial expansion).

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_lascoux"></a>

#### transition\_lascoux

```python
@classmethod
def transition_lascoux(cls, key)
```

Expand a word in `LascouxBasis`: the K-theoretic analogue of `transition_key` using WC graphs.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_glide"></a>

#### transition\_glide

```python
@classmethod
def transition_glide(cls, key)
```

Expand a word in `GlideBasis`: for each WC graph of weight ``key``, take the length vector of
its ``dst`` (destandardization); the first graph seen at each weight is the representative
and only graphs with that same ``dst`` contribute.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
@classmethod
def transition_fundamental_slide(cls, key)
```

Expand a word in `FundamentalSlideBasis` as the transpose of the polynomial side: the
coefficient of ``candidate`` is the coefficient of the monomial ``x^key`` in the fundamental
slide polynomial of ``candidate``, over all weak compositions of the same length and size.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, key)
```

Transition a word key (composition) to the Grothendieck basis.

Coefficient of ``G_w`` is the number of WC graphs of permutation ``w``
and weight ``key``, multiplied by ``beta^(|key|-inv(w))``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition"></a>

#### transition

```python
@classmethod
@cache
def transition(cls, other_basis)
```

Key -> ``{key: coeff}`` function into ``other_basis``; dispatches to the ``transition_*``
method for each directly supported basis and otherwise goes through `SchubertBasis`.

<a id="schubmult.rings.free_algebra.z_basis"></a>

# schubmult.rings.free\_algebra.z\_basis

`ZBasis`: free-algebra basis indexed by compositions with no zeros, related to `SchubertBasis`
by shifting code entries by one and dropping zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis"></a>

## ZBasis Objects

```python
class ZBasis(FreeAlgebraBasis)
```

Z basis of the free algebra.

Keys are tuples of positive integers (no zeros). The Z basis is
related to the Schubert basis by incrementing/decrementing code
entries by 1 and dropping zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a Z basis key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.pare_schubert"></a>

#### pare\_schubert

```python
@staticmethod
def pare_schubert(perm)
```

Extract the nonzero trimcode of *perm*, or None if it contains interior zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two Z basis keys via shifted Schubert multiplication.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Z``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ZBasis to *other_basis*.

<a id="schubmult.rings.nsym"></a>

# schubmult.rings.nsym

`NSym`: a `FreeAlgebra` indexed by compositions, multiplied through the separated-descents Schubert ring.

A key is a composition ``alpha`` (a tuple of positive integers), printed ``N(alpha)``, and is
identified with the Schubert key ``(uncode(alpha - 1), len(alpha))`` of
`schubmult.rings.schubert.separated_descents.SeparatedDescentsRing` via `NSym.sepify` /
`NSym.from_sep`. The product is the separated-descents product transported back to
compositions; when ``FreeAlgebra.CAP`` is set the result is truncated to keys of at most that
length. Right multiplication by a Schubert element acts by the skew operation ``/``.

<a id="schubmult.rings.nsym.NSym"></a>

## NSym Objects

```python
class NSym(FreeAlgebra)
```

Free algebra on compositions with the separated-descents product. See the module docstring.

<a id="schubmult.rings.nsym.NSym.__init__"></a>

#### \_\_init\_\_

```python
def __init__(domain=None)
```

Create the ring over ``domain`` (default ``EXRAW``); the empty composition is the identity.

<a id="schubmult.rings.nsym.NSym.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Display as ``N(alpha)``.

<a id="schubmult.rings.nsym.NSym.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Scale coefficients by the scalar ``other``.

<a id="schubmult.rings.nsym.NSym.sepify"></a>

#### sepify

```python
def sepify(elem)
```

Map ``alpha -> (uncode(alpha - 1), len(alpha))`` into the separated-descents Schubert ring.

<a id="schubmult.rings.nsym.NSym.from_sep"></a>

#### from\_sep

```python
def from_sep(elem)
```

Inverse of `sepify`: pad or cut the code of ``perm`` to length ``n`` and add 1 to each entry.

<a id="schubmult.rings.nsym.NSym.mul"></a>

#### mul

```python
def mul(elem, other)
```

Scalar multiplication, or the separated-descents product of two elements (truncated by
``FreeAlgebra.CAP`` if set).

<a id="schubmult.rings.nsym.NSymElement"></a>

## NSymElement Objects

```python
class NSymElement(FreeAlgebraElement)
```

Element of `NSym`: a dict from compositions to coefficients with SymPy-compatible printing.

<a id="schubmult.rings.nsym.NSymElement.__rmul__"></a>

#### \_\_rmul\_\_

```python
def __rmul__(other)
```

Scalar on the left, or a Schubert element acting by the skew operation ``self / perm``.

<a id="schubmult.rings.polynomial_algebra"></a>

# schubmult.rings.polynomial\_algebra

Polynomial algebra ring and pre-built basis instances.

Exports the core :class:`PolynomialAlgebra` and :class:`PolynomialAlgebraElement`
classes as well as ready-to-use ring instances:

- ``Schub`` — Schubert polynomial basis
- ``Forest`` — Forest polynomial basis
- ``Key`` — Key polynomial (Demazure character) basis
- ``FSlide`` — Fundamental slide polynomial basis
- ``Monomial`` — Standard monomial basis

<a id="schubmult.rings.polynomial_algebra._core"></a>

# schubmult.rings.polynomial\_algebra.\_core

`PolynomialAlgebra`: the polynomial ring ``Z[x_1, x_2, ...]`` with a pluggable basis.

The ring itself is basis-agnostic; a `PolynomialBasis` instance supplies the key
type, the product rule, the coproduct, and the transitions to/from the monomial
basis. Elements of rings with different bases are interconverted via
``change_basis``. ``PA`` is the standard monomial-basis instance in ``x``; the
pre-built instances for other bases (``Schub``, ``Key``, ``FSlide``, ...) live in
the package ``__init__``.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement"></a>

## PolynomialAlgebraElement Objects

```python
class PolynomialAlgebraElement(BaseRingElement)
```

Element of a polynomial algebra, stored as a dict mapping basis keys to coefficients.

Keys are exponent tuples (in the monomial basis) or basis-specific keys
depending on the parent ring's basis. Supports arithmetic, basis changes,
and duality pairing with free algebra elements.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

Return a dict mapping printing terms to sympified coefficients.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.branch"></a>

#### branch

```python
def branch(index)
```

Split the variables at ``index``: ``x_1..x_index`` on the left tensor factor, the rest on the
right, returned in the tensor square of this ring's basis.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Sum of ``branch(index)`` over every split point (the full variable-splitting coproduct).

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis: type)
```

Convert this element to another polynomial basis.

**Arguments**:

- `other_basis` - A basis class, basis instance, or callable returning a basis.
  

**Returns**:

  A new PolynomialAlgebraElement in the target basis's ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.expand"></a>

#### expand

```python
def expand()
```

Expand this element into an explicit polynomial expression.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.apply_dual_element"></a>

#### apply\_dual\_element

```python
def apply_dual_element(dual_elem)
```

Pair this polynomial element with a dual free algebra element.

Converts *self* to the monomial basis and *dual_elem* to the word
basis, then sums products of matching coefficients.

**Arguments**:

- `dual_elem` - A FreeAlgebraElement to pair with.
  

**Returns**:

  The scalar pairing value.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra"></a>

## PolynomialAlgebra Objects

```python
class PolynomialAlgebra(BaseRing)
```

Polynomial algebra ring with a configurable basis.

The algebra operates on :class:`PolynomialAlgebraElement` instances whose
keys are determined by the chosen basis. Supports multiplication, basis
changes, coproducts, and conversion from symbolic expressions.

**Arguments**:

- `basis` - A basis instance (e.g. ``MonomialBasis(x)``).
- `domain` - Coefficient domain (default ``EXRAW``).

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.__init__"></a>

#### \_\_init\_\_

```python
def __init__(basis, domain=None)
```

Initialize a PolynomialAlgebra with the given basis and coefficient domain.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.genset"></a>

#### genset

```python
@property
def genset()
```

The basis's generating set.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(key)
```

Compute the coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via the basis product rule.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.new"></a>

#### new

```python
def new(*x)
```

Create a new element from the given key or expression.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.from_expr"></a>

#### from\_expr

```python
def from_expr(x, length=None)
```

Create an element from a symbolic expression.

Parses *x* into monomials, then transitions to this ring's basis.

**Arguments**:

- `x` - A symbolic polynomial expression.
- `length` - Optional fixed number of variables.
  

**Returns**:

  A PolynomialAlgebraElement in this ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Return the display symbol for basis key *k*.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.from_dict"></a>

#### from\_dict

```python
def from_dict(element)
```

Construct an element from a dict of ``{key: coefficient}`` pairs.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce a raw value into the coefficient domain.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.anti\_schubert\_poly\_basis

`AntiSchubertPolyBasis`: the anti-Schubert (``w0``-conjugated Schubert) polynomial basis of `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis"></a>

## AntiSchubertPolyBasis Objects

```python
class AntiSchubertPolyBasis(PolynomialBasis)
```

Anti-Schubert polynomial basis.

Keys are ``(Permutation, length)`` pairs. This basis reverses the
monomial ordering relative to the standard Schubert basis, with
the coproduct correspondingly reversed.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the reversed coproduct of an anti-Schubert key.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two anti-Schubert keys using the underlying Schubert ring.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_key_key"></a>

#### transition\_key\_key

```python
def transition_key_key(key)
```

Decompose an anti-Schubert polynomial into key polynomials with reversed weights.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_key"></a>

#### transition\_key

```python
def transition_key(dct)
```

Transition an anti-Schubert dict to the key polynomial basis.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand an anti-Schubert key into reversed monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_forest"></a>

#### transition\_forest

```python
def transition_forest(dct)
```

Transition an anti-Schubert dict to the forest polynomial basis.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from anti-Schubert basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis"></a>

# schubmult.rings.polynomial\_algebra.base\_polynomial\_basis

`PolynomialBasis`: the abstract interface a basis must implement to plug into `PolynomialAlgebra`.

A basis defines its key type (``is_key``/``as_key``/``zero_monom``), how to print a
key, and ``transition(other_basis)`` -- a function converting coefficient dicts
into another basis. Products, coproducts, expansion, and parsing from expressions
all have default implementations that route through the `MonomialBasis`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis"></a>

## PolynomialBasis Objects

```python
class PolynomialBasis(ABC)
```

Abstract base class for polynomial algebra bases.

Subclasses define how keys are represented, how to transition between
bases, and how to expand elements into explicit polynomials. Default
implementations delegate through the :class:`MonomialBasis`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.genset"></a>

#### genset

```python
@property
def genset()
```

The generating set (variable alphabet).

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.is_key"></a>

#### is\_key

```python
@abstractmethod
def is_key(x)
```

Return True if *x* is a valid key for this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.as_key"></a>

#### as\_key

```python
@abstractmethod
def as_key(x)
```

Normalize *x* into a canonical key for this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.attach_key"></a>

#### attach\_key

```python
def attach_key(dct)
```

Normalize all keys in *dct* via :meth:`as_key`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.zero_monom"></a>

#### zero\_monom

```python
@property
@abstractmethod
def zero_monom()
```

The key of the multiplicative identity.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.monomial_basis"></a>

#### monomial\_basis

```python
@property
def monomial_basis()
```

The `MonomialBasis` over the same generating set (the hub for default transitions).

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.transition"></a>

#### transition

```python
@abstractmethod
def transition(other_basis)
```

Return a function mapping dicts of this basis to dicts in *other_basis*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.from_expr"></a>

#### from\_expr

```python
def from_expr(expr, length=None)
```

Parse a symbolic expression into this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.printing_term"></a>

#### printing\_term

```python
@abstractmethod
def printing_term(k)
```

Return the display symbol for key *k*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.compose_transition"></a>

#### compose\_transition

```python
@staticmethod
def compose_transition(tkeyfunc, output)
```

Apply a transition function to a dict of basis elements.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.change_tensor_basis"></a>

#### change\_tensor\_basis

```python
@classmethod
def change_tensor_basis(cls, tensor_elem, basis1, basis2)
```

Change the bases of both factors of a tensor element.

**Arguments**:

- `tensor_elem` - An element of a tensor product ring.
- `basis1` - Target basis for the left factor.
- `basis2` - Target basis for the right factor.
  

**Returns**:

  The tensor element re-expressed in the new bases.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a basis dict into an explicit polynomial expression.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the coproduct of *key* by delegating through the monomial basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to the monomial basis and back.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class.

<a id="schubmult.rings.polynomial_algebra.composition_schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.composition\_schubert\_poly\_basis

`CompositionSchubertPolyBasis`: Schubert polynomials re-indexed by weak compositions (Lehmer codes)
instead of permutations, wrapping `SchubertPolyBasis`.

<a id="schubmult.rings.polynomial_algebra.composition_schubert_poly_basis.CompositionSchubertPolyBasis"></a>

## CompositionSchubertPolyBasis Objects

```python
class CompositionSchubertPolyBasis(SchubertPolyBasis)
```

Wrapper basis for Schubert polynomials indexed by weak compositions.

Keys are weak compositions interpreted as Lehmer codes. The underlying
computations are delegated to :class:`SchubertPolyBasis`, while display
uses the same Schubert printing as the corresponding permutation key.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.double\_forest\_poly\_basis

`DoubleForestPolyBasis`: the double (two-alphabet) forest polynomial basis of `PolynomialAlgebra`;
see `schubmult.combinatorics.double_forest` for the underlying polynomials.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis"></a>

## DoubleForestPolyBasis Objects

```python
class DoubleForestPolyBasis(PolynomialBasis)
```

Abstract double forest polynomial basis.

Keys are weak compositions indexing double forest basis elements DF[key],
with polynomial coefficients in a second generating set (equivariant vars).

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis.basis_polynomial"></a>

#### basis\_polynomial

```python
def basis_polynomial(key)
```

Expand one double-forest basis key as a polynomial in x,t.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis.basis_forest_expansion"></a>

#### basis\_forest\_expansion

```python
def basis_forest_expansion(key, length)
```

Expand one double-forest basis key into ForestPolyBasis in x.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.elem\_sym\_poly\_basis

`ElemSymPolyBasis`: the basis of products of elementary symmetric polynomials ``e_p(x_1..x_k)`` for `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis"></a>

## ElemSymPolyBasis Objects

```python
class ElemSymPolyBasis(PolynomialBasis)
```

Elementary symmetric polynomial basis.

Keys are tuples encoding products of elementary symmetric polynomials
e_k(x_1, ..., x_n). Each key specifies degrees and variable counts
for the elementary symmetric factors.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition from elementary symmetric basis to Schubert basis.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from elementary symmetric basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from this basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.forest\_poly\_basis

`ForestPolyBasis`: the forest polynomial basis (Nadeau-Spink-Tewari) of `PolynomialAlgebra`, indexed
by weak compositions via indexed forests (`schubmult.combinatorics.indexed_forests`).

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis"></a>

## ForestPolyBasis Objects

```python
class ForestPolyBasis(PolynomialBasis)
```

Forest polynomial basis.

Keys are weak compositions encoding indexed forests. Forest polynomials
are computed by summing over decreasing labelings of the corresponding
forest structure.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a forest key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a forest basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from forest basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
def transition_fundamental_slide(dct)
```

Transition from forest basis to fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.to_fundamental_slide"></a>

#### to\_fundamental\_slide

```python
def to_fundamental_slide(key)
```

Express a single forest key in the fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`ForestBasis`).

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from forest basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two forest keys by transitioning through the Schubert basis.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.fundamental\_slide\_poly\_basis

`FundamentalSlidePolyBasis`: the fundamental slide polynomial basis (Assaf-Searles) of `PolynomialAlgebra`, indexed by weak compositions.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.get_descent_composition"></a>

#### get\_descent\_composition

```python
def get_descent_composition(word)
```

Compute the descent composition of a word.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.slide_product"></a>

#### slide\_product

```python
def slide_product(a, b)
```

Compute the structure constants for multiplying two fundamental slide polynomials.

Given weak compositions *a* and *b*, returns a dict mapping result
compositions to their coefficients in the fundamental slide expansion
of the product.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis"></a>

## FundamentalSlidePolyBasis Objects

```python
class FundamentalSlidePolyBasis(PolynomialBasis)
```

Fundamental slide polynomial basis.

Keys are weak compositions. Fundamental slide polynomials provide a
basis that refines Schubert polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a slide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`FundamentalSlideBasis`).

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a slide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from fundamental slide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from fundamental slide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two fundamental slide keys using the slide product rule.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.glide\_poly\_basis

`GlidePolyBasis`: the glide polynomial basis (K-theoretic analogue of fundamental slides) of `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.glide_monomials"></a>

#### glide\_monomials

```python
def glide_monomials(key)
```

Monomial expansion of the glide polynomial :math:`\mathcal{G}_{key}`.

Returns a dict mapping each exponent tuple ``v`` (a weak composition of the
same length as ``key``) to the integer coefficient of the monomial
:math:`x^v`, i.e. the number of glides of ``key`` with weight ``v``. The
corresponding power of ``beta`` for the weight ``v`` is
``sum(v) - sum(key)`` (the excess), which is constant across all glides of a
given weight, so it need not be stored explicitly.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.glide_product"></a>

#### glide\_product

```python
def glide_product(key1, key2)
```

Structure constants for a product of two glide polynomials.

Implements the Littlewood-Richardson rule of O. Pechenik and D. Searles,
"Decompositions of Grothendieck Polynomials" (arXiv:1611.02545), Theorem
4.9, which expands the product of the glide polynomials indexed by the weak
compositions ``key1`` and ``key2`` in the glide basis:

.. math::

    \mathcal{G}_a \, \mathcal{G}_b
        = \sum_c \beta^{|c| - |a| - |b|} \, g_{a,b}^{c} \, \mathcal{G}_c .

Rather than enumerating the genomic shuffle set directly, we compute the
(uniquely determined) coefficients by expanding the product in monomials and
straightening into the glide basis with the leading-term algorithm from the
proof that the glide polynomials form a basis (Theorem 2.6). Because the
excess of a glide of ``v`` equals ``sum(v) - sum(index)``, the power of
``beta`` is recovered from the total degree and only the positive integer
multiplicities :math:`g_{a,b}^{c}` are returned.

Both compositions are padded with trailing zeros to a common length ``n``;
every key ``c`` in the returned dict is a weak composition of length ``n``.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis"></a>

## GlidePolyBasis Objects

```python
class GlidePolyBasis(PolynomialBasis)
```

Glide polynomial basis.

Keys are weak compositions. Glide polynomials provide a
basis that refines Grothendieck polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a glide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GlideBasis`).

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a glide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from glide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from glide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two glide keys using the glide product rule.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.grothendieck\_poly\_basis

`GrothendieckPolyBasis`: the Grothendieck polynomial basis of `PolynomialAlgebra`, at the
specialization ``beta = 1`` (without loss of generality: ``beta`` is recovered from the grading,
since the degree ``inv(w) + d`` part of ``G_w`` carries ``beta^d``).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis"></a>

## GrothendieckPolyBasis Objects

```python
class GrothendieckPolyBasis(PolynomialBasis)
```

Grothendieck polynomial basis at ``beta = 1``.

Keys are ``(Permutation, length)`` pairs. Grothendieck polynomials form
the canonical basis for the polynomial algebra in Grothendieck calculus,
dual to the :class:`GrothendieckBasis` of the free algebra. The ``beta`` parameter
is set to 1 without loss of generality (see the module docstring).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two Grothendieck keys using the Grothendieck ring multiplication.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition from Grothendieck basis to separated descents basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_glide_key"></a>

#### transition\_glide\_key

```python
def transition_glide_key(key)
```

Decompose a Grothendieck polynomial into glide polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_glide"></a>

#### transition\_glide

```python
def transition_glide(dct)
```

Transition a Grothendieck dict to the glide polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a Grothendieck key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GrothendieckBasis`).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_grove_key"></a>

#### transition\_grove\_key

```python
def transition_grove_key(key)
```

Decompose a Grothendieck polynomial into grove polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_grove"></a>

#### transition\_grove

```python
def transition_grove(dct)
```

Transition a Grothendieck dict to the grove polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_lascoux_key"></a>

#### transition\_lascoux\_key

```python
def transition_lascoux_key(key)
```

Decompose a Grothendieck polynomial into Lascoux polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_lascoux"></a>

#### transition\_lascoux

```python
def transition_lascoux(dct)
```

Transition a Grothendieck dict to the Lascoux polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Grothendieck basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.grove\_poly\_basis

`GrovePolyBasis`: the grove polynomial basis (K-theoretic analogue of forest polynomials) of `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis"></a>

## GrovePolyBasis Objects

```python
class GrovePolyBasis(PolynomialBasis)
```

Grove polynomial basis.

Keys are weak compositions encoding indexed forests. Grove polynomials are
the ``beta``-deformed (set-valued) forest polynomials, computed by summing
over set-valued labelings of the corresponding forest structure.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a grove key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a grove basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from grove basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two grove keys using the glide product rule.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from grove basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GroveBasis`).

<a id="schubmult.rings.polynomial_algebra.key_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.key\_poly\_basis

`KeyPolyBasis`: the key polynomial (Demazure character) basis of `PolynomialAlgebra`, indexed by weak compositions.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.traverse_demaz"></a>

#### traverse\_demaz

```python
def traverse_demaz(pl, w)
```

Traverse the Demazure graph starting from plactic element *pl* along the code word of *w*.

Yields all distinct elements reachable by successive lowering operators.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis"></a>

## KeyPolyBasis Objects

```python
class KeyPolyBasis(PolynomialBasis)
```

Key polynomial (Demazure character) basis.

Keys are weak compositions. Key polynomials are characters of Demazure
modules, computed by applying Demazure operators to highest-weight
plactic tableaux.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`KeyBasis`).

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a key polynomial key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a key basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from key basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from key basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.lascoux\_poly\_basis

`LascouxPolyBasis`: the Lascoux polynomial (K-theoretic key polynomial) basis of `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis"></a>

## LascouxPolyBasis Objects

```python
class LascouxPolyBasis(PolynomialBasis)
```

Lascoux polynomial basis.

Keys are weak compositions. Lascoux polynomials provide a
basis that refines Grothendieck polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a glide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`LascouxBasis`).

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a Lascoux basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from Lascoux basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_glide_key"></a>

#### transition\_glide\_key

```python
def transition_glide_key(key)
```

Transition a Lascoux key to the glide basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_glide"></a>

#### transition\_glide

```python
def transition_glide(dct)
```

Transition from Lascoux basis to glide basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Lascoux basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two Lascoux keys using the Lascoux product rule.

<a id="schubmult.rings.polynomial_algebra.monomial_basis"></a>

# schubmult.rings.polynomial\_algebra.monomial\_basis

`MonomialBasis`: the standard monomial basis ``x^a`` (keys are exponent tuples) of `PolynomialAlgebra`.

This is the hub basis: every other `PolynomialBasis` transitions through it by default.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis"></a>

## MonomialBasis Objects

```python
class MonomialBasis(PolynomialBasis)
```

Standard monomial basis for the polynomial algebra.

Keys are tuples of nonnegative integers representing exponent vectors.
This is the fundamental basis through which other bases transition
by default, and is dual to the :class:`WordBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the deconcatenation coproduct on a monomial key.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two monomial keys by component-wise addition of exponents.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.expand_monom"></a>

#### expand\_monom

```python
def expand_monom(monom)
```

Convert an exponent tuple to a monomial expression in the generating set.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a dict of monomial keys into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_slide"></a>

#### transition\_slide

```python
def transition_slide(dct, other_basis)
```

Transition a monomial dict to a slide-type basis via triangular inversion.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_slide_monom"></a>

#### transition\_slide\_monom

```python
def transition_slide_monom(other_basis, monom, coeff=S.One)
```

Express a single monomial in a slide-type basis by dominance-order inversion.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_anti_schubert"></a>

#### transition\_anti\_schubert

```python
def transition_anti_schubert(dct, other_basis)
```

Transition monomials to the anti-Schubert basis by reversing exponents.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition monomials to the Schubert basis by grouping by length.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_double_forest"></a>

#### transition\_double\_forest

```python
def transition_double_forest(dct, other_basis)
```

Transition monomials to DoubleForestPolyBasis via ForestPolyBasis.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`WordBasis`).

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from monomial basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.from_expr"></a>

#### from\_expr

```python
def from_expr(expr, length=None)
```

Parse a symbolic expression into monomial-basis coefficient dict.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.monomial\_slide\_poly\_basis

`MonomialSlidePolyBasis`: the monomial slide polynomial basis (Assaf-Searles) of `PolynomialAlgebra`.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis"></a>

## MonomialSlidePolyBasis Objects

```python
class MonomialSlidePolyBasis(PolynomialBasis)
```

Monomial slide polynomial basis.

Keys are weak compositions. Monomial slide polynomials refine
key polynomials and are coarser than monomials, using a recursive
construction based on the first nonzero entry.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a monomial slide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from monomial slide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a monomial slide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from monomial slide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.polynomial_basis"></a>

# schubmult.rings.polynomial\_algebra.polynomial\_basis

Re-export hub for all polynomial basis classes.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.schubert\_poly\_basis

`SchubertPolyBasis`: the Schubert polynomial basis of `PolynomialAlgebra`, indexed by permutations.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis"></a>

## SchubertPolyBasis Objects

```python
class SchubertPolyBasis(PolynomialBasis)
```

Schubert polynomial basis.

Keys are ``(Permutation, length)`` pairs. Schubert polynomials form
the canonical basis for the polynomial algebra in Schubert calculus,
dual to the :class:`SchubertBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the coproduct of a Schubert key by splitting variable sets.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two Schubert keys using the Schubert ring multiplication.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
def transition_grothendieck(dct)
```

Transition a Schubert dict to the Grothendieck polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_sepdesc"></a>

#### transition\_sepdesc

```python
def transition_sepdesc(dct, other_basis)
```

Transition from Schubert basis to separated descents basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_elementary"></a>

#### transition\_elementary

```python
def transition_elementary(dct, other_basis)
```

Transition from Schubert basis to elementary symmetric basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key_fundamental_slide"></a>

#### transition\_key\_fundamental\_slide

```python
def transition_key_fundamental_slide(perm, n)
```

Decompose a Schubert polynomial into fundamental slide polynomials via quasi-Yamanouchi RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
def transition_fundamental_slide(dct)
```

Transition a Schubert dict to the fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key_key"></a>

#### transition\_key\_key

```python
def transition_key_key(key)
```

Decompose a Schubert polynomial into key polynomials via highest-weight RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key"></a>

#### transition\_key

```python
def transition_key(dct)
```

Transition a Schubert dict to the key polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a Schubert key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`SchubertBasis`).

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_forest_key"></a>

#### transition\_forest\_key

```python
def transition_forest_key(key)
```

Decompose a Schubert polynomial into forest polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_forest"></a>

#### transition\_forest

```python
def transition_forest(dct)
```

Transition a Schubert dict to the forest polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Schubert basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.sepdesc\_poly\_basis

`SepDescPolyBasis`: the separated-descents polynomial basis of `PolynomialAlgebra`, indexed by
``(perm, num_vars)`` pairs (see `schubmult.rings.schubert.separated_descents`).

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis"></a>

## SepDescPolyBasis Objects

```python
class SepDescPolyBasis(PolynomialBasis)
```

Separated descents polynomial basis.

Keys are ``(Permutation, Permutation, k)`` triples. Elements are
products of pairs of Schubert polynomials, parameterized by a
separation level *k*.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two separated-descents keys by transitioning through Schubert.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct, other_basis)
```

Transition from separated descents to Schubert basis by multiplying the pair factors.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from separated descents to *other_basis*.

<a id="schubmult.rings.printing"></a>

# schubmult.rings.printing

SymPy atoms used to display ring basis elements.

Every ring's ``printing_term(key)`` returns a `PrintingTerm` subclass instance: an inert SymPy
``Expr`` atom (``args == ()``, so SymPy never traverses into it) that knows how to render itself
for ``str``, pretty printing, and LaTeX. Instances are interned via cached ``__xnew_cached__``
constructors so equal keys give identical objects. The subclasses cover single/double Schubert
(``S``/``DS``), quantum (``QS``/``QDS``, ``QPS``/``QPDS``), Grothendieck (``G``/``DG``), separated
descents (``Xi``), and a `GenericPrintingTerm` ``name(key)`` fallback.

<a id="schubmult.rings.printing.PrintingTerm"></a>

## PrintingTerm Objects

```python
class PrintingTerm(ssymb.Expr)
```

Base display atom carrying a key, generating set, coefficient generating set, and prefix.

<a id="schubmult.rings.printing.GenericPrintingTerm"></a>

## GenericPrintingTerm Objects

```python
class GenericPrintingTerm(PrintingTerm)
```

Displays a key as ``name(key)`` (e.g. ``AGx(perm, n)``, ``N(2, 1)``); the identity key prints as ``1``.

<a id="schubmult.rings.printing.TypedPrintingTerm"></a>

## TypedPrintingTerm Objects

```python
class TypedPrintingTerm(PrintingTerm)
```

Displays a key by delegating to the key's own printer (used for keys that are themselves
printable objects such as RC graphs).

<a id="schubmult.rings.printing.DSchubPoly"></a>

## DSchubPoly Objects

```python
class DSchubPoly(PrintingTerm)
```

Schubert polynomial term: ``S<genset>(perm)`` or ``DS<genset>(perm, <coeff_genset>)``.

<a id="schubmult.rings.printing.SepDescSchubPoly"></a>

## SepDescSchubPoly Objects

```python
class SepDescSchubPoly(PrintingTerm)
```

Separated-descents term for the key ``(perm, numvars)``: ``Xi_{perm}^{numvars}``.

<a id="schubmult.rings.printing.QDSchubPoly"></a>

## QDSchubPoly Objects

```python
class QDSchubPoly(PrintingTerm)
```

Quantum Schubert term: ``QS<genset>(perm)`` or ``QDS<genset>(perm, <coeff_genset>)``.

<a id="schubmult.rings.printing.PQDSchubPoly"></a>

## PQDSchubPoly Objects

```python
class PQDSchubPoly(PrintingTerm)
```

Parabolic quantum Schubert term, tagged with the index composition:
``QPS<genset>(comp)(perm)`` or ``QPDS<genset>(comp)(perm, <coeff_genset>)``.

<a id="schubmult.rings.printing.GrothendieckPoly"></a>

## GrothendieckPoly Objects

```python
class GrothendieckPoly(PrintingTerm)
```

Grothendieck term ``G<genset>(perm)``, or ``G<genset>(perm, numvars)`` when the key carries a variable count.

<a id="schubmult.rings.printing.DoubleGrothendieckPoly"></a>

## DoubleGrothendieckPoly Objects

```python
class DoubleGrothendieckPoly(PrintingTerm)
```

Double Grothendieck term ``DG<genset>(perm, <coeff_genset>)``.

<a id="schubmult.rings.product_ring"></a>

# schubmult.rings.product\_ring

`ProductRing`: array-backed componentwise product of Schubert-family rings (draft).

An element holds a numpy object array ``_arr`` with one ring element per factor; arithmetic is
componentwise on the array. This is an older sketch of the same idea as
`schubmult.rings.direct_product_ring.DirectProductRing`, which is the supported implementation;
`ProductRing` is not exported from the package.

<a id="schubmult.rings.product_ring.ProductRing"></a>

## ProductRing Objects

```python
class ProductRing(BaseSchubertRing)
```

Componentwise product of rings backed by a numpy array of factor elements. See the module docstring.

<a id="schubmult.rings.product_ring.ProductRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(*rings)
```

Flatten nested product rings and pool the factors' generators and coefficient generators.

<a id="schubmult.rings.product_ring.ProductRing.rings"></a>

#### rings

```python
@property
def rings()
```

The (flattened) tuple of factor rings.

<a id="schubmult.rings.product_ring.ProductRing.new"></a>

#### new

```python
def new(x)
```

Wrap a sequence of factor elements as an element of this ring.

<a id="schubmult.rings.product_ring.ProductRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(*x)
```

Build an element from one input per factor (or a single sequence/element), coercing each
input in its factor ring.

<a id="schubmult.rings.product_ring.ProductBasisElement"></a>

## ProductBasisElement Objects

```python
class ProductBasisElement(PrintingTerm)
```

Printing term for a `ProductRing` element; renders the factors joined by ``#``.

<a id="schubmult.rings.product_ring.ProductRingElement"></a>

## ProductRingElement Objects

```python
class ProductRingElement(BaseSchubertElement)
```

Element of a `ProductRing`; arithmetic operators act componentwise on ``_arr``.

<a id="schubmult.rings.quasisymmetric_functions"></a>

# schubmult.rings.quasisymmetric\_functions

`QSym`: quasisymmetric functions in the monomial basis ``M_alpha``.

Keys are compositions; the product is the quasi-shuffle (stuffle) of compositions, and
``expand(n)`` gives the monomial quasisymmetric polynomial in ``n`` variables. `QSym.quasi_schur`
builds quasi-Schur functions by enumerating standard composition tableaux.

<a id="schubmult.rings.quasisymmetric_functions.monomial_quasisym"></a>

#### monomial\_quasisym

```python
def monomial_quasisym(comp, length, genset)
```

The monomial quasisymmetric polynomial ``M_comp(x_1, ..., x_length)``: the sum of
``x_{i_1}^{c_1} ... x_{i_k}^{c_k}`` over ``i_1 < ... < i_k <= length``, built by recursion on
whether ``x_length`` is used. Zero if ``comp`` contains a zero part.

<a id="schubmult.rings.quasisymmetric_functions.stuffle"></a>

#### stuffle

```python
def stuffle(alpha, beta)
```

The quasi-shuffle (stuffle) product of two compositions: at each step take the first part
of ``alpha``, the first part of ``beta``, or their sum. Returns ``{composition: coeff}``;
this is the product rule of the monomial basis ``M_alpha M_beta``.

<a id="schubmult.rings.quasisymmetric_functions.quasi_schur_to_monomial"></a>

#### quasi\_schur\_to\_monomial

```python
def quasi_schur_to_monomial(comp)
```

Monomial-basis expansion of the quasi-Schur function of shape ``comp``: counts standard
composition tableaux (rows strictly increasing, columns weakly increasing) of that shape by
the descent composition of their row reading word. Enumerates all ``n!`` fillings, so only
small shapes are practical.

<a id="schubmult.rings.quasisymmetric_functions.QSymElement"></a>

## QSymElement Objects

```python
class QSymElement(BaseSchubertElement)
```

Element of `QSym`: a dict from compositions to coefficients in the monomial basis.

<a id="schubmult.rings.quasisymmetric_functions.QSymElement.expand"></a>

#### expand

```python
def expand(num_vars)
```

The quasisymmetric polynomial in ``num_vars`` variables of the ring's generating set.

<a id="schubmult.rings.quasisymmetric_functions.QSym"></a>

## QSym Objects

```python
class QSym(BaseSchubertRing)
```

Quasisymmetric functions in the monomial basis; ``QSym()(2, 1)`` is ``M_(2,1)``. See the module docstring.

<a id="schubmult.rings.quasisymmetric_functions.QSym.mul_pair"></a>

#### mul\_pair

```python
def mul_pair(a, b)
```

Product of two basis compositions: the `stuffle`.

<a id="schubmult.rings.quasisymmetric_functions.QSym.mul"></a>

#### mul

```python
def mul(a, b)
```

Bilinear extension of `mul_pair`.

<a id="schubmult.rings.quasisymmetric_functions.QSym.printing_term"></a>

#### printing\_term

```python
def printing_term(comp)
```

Display as ``Mx(alpha)`` (label from the generating set).

<a id="schubmult.rings.quasisymmetric_functions.QSym.new"></a>

#### new

```python
def new(*x)
```

The basis element ``M_x`` for the composition given as positional parts.

<a id="schubmult.rings.quasisymmetric_functions.QSym.quasi_schur"></a>

#### quasi\_schur

```python
def quasi_schur(*comp)
```

The quasi-Schur function of shape ``comp`` in the monomial basis (see `quasi_schur_to_monomial`).

>>> QS = QSym()
>>> QS.quasi_schur(2, 1)

<a id="schubmult.rings.schubert"></a>

# schubmult.rings.schubert

Schubert-family rings: the main user-facing algebra interface.

The most commonly used entry points are the ring *instances*:

- ``Sx``: ordinary (single) Schubert polynomials, ``Sx([3, 1, 2]) * Sx([2, 1, 3])``.
- ``DSx``: double Schubert polynomials (second alphabet ``y``).
- ``Gx`` / ``DGx``: (double) Grothendieck polynomials.
- ``QSx`` / ``QDSx``: quantum (double) Schubert polynomials.
- ``QPSx`` / ``QPDSx``: parabolic quantum (double) Schubert polynomials.

Each instance is an object of the corresponding ``*Ring`` class; calling it with
a permutation (or Lehmer code, or a polynomial expression) yields a ``*Element``
that supports ``+``, ``*``, ``.expand()``, and conversion between bases. All ring
classes derive from `BaseSchubertRing` and dispatch their products to the kernels
in `schubmult.mult`.

Everything here is imported lazily (see ``__getattr__``) so ``import schubmult``
stays fast.

<a id="schubmult.rings.schubert.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

Lazily import and cache the requested export from its defining submodule.

<a id="schubmult.rings.schubert.__dir__"></a>

#### \_\_dir\_\_

```python
def __dir__()
```

Include lazily-exported names in ``dir()``.

<a id="schubmult.rings.schubert.base_schubert_ring"></a>

# schubmult.rings.schubert.base\_schubert\_ring

Abstract base classes shared by every Schubert-family ring.

`BaseSchubertRing` holds a primary generating set (``genset``, the ``x`` variables)
and a coefficient generating set (``coeff_genset``, the ``y``/``z`` variables, or
``None`` for single Schubert polynomials), and defines the ring-level hooks that
concrete rings fill in: which multiplication kernel to use, how to expand a basis
element to a polynomial, how to print it, and how to change basis. `BaseSchubertElement`
is the corresponding dict-like element type (``{Permutation: coefficient}``).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement"></a>

## BaseSchubertElement Objects

```python
class BaseSchubertElement(BaseRingElement)
```

A linear combination of Schubert-family basis elements, stored as ``{Permutation: coeff}``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.mult_poly"></a>

#### mult\_poly

```python
def mult_poly(poly)
```

Multiply this element by an arbitrary polynomial ``poly`` in the ring's ``genset`` variables.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.in_schubert_schur_basis"></a>

#### in\_schubert\_schur\_basis

```python
def in_schubert_schur_basis(numvars)
```

Expand into the Schubert-tensor-Schur basis of the tensor square ring, splitting off the
symmetric part in the last ``numvars`` variables.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.in_SEM_basis"></a>

#### in\_SEM\_basis

```python
def in_SEM_basis(elem_func=None)
```

Expand as a polynomial in elementary symmetric functions (the "SEM" presentation), using
``elem_func`` (default: the ring's symbolic ``symbol_elem_func``) as the elementary symmetric symbol.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms ``coeff * basis_symbol`` sorted by permutation length then lexicographically (sympy printing hook).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

With ``deep=True`` (default) expand to an explicit polynomial in the variables; with
``deep=False`` only expand each coefficient, keeping the Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_expr"></a>

#### as\_expr

```python
def as_expr()
```

Sum of the ``as_terms()`` as a sympy ``Add``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Expand to an explicit polynomial: ``sum coeff * SchubertPoly(perm)``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_classical"></a>

#### as\_classical

```python
def as_classical()
```

Re-express in the classical (non-quantum) Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.as_quantum"></a>

#### as\_quantum

```python
def as_quantum()
```

Re-express in the quantum Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.almosteq"></a>

#### almosteq

```python
def almosteq(other)
```

Equality up to coefficient expansion (handles elements of different but compatible rings).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertElement.strip_zeros"></a>

#### strip\_zeros

```python
def strip_zeros()
```

Drop basis elements whose coefficient is exactly zero.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing"></a>

## BaseSchubertRing Objects

```python
class BaseSchubertRing(BaseRing)
```

Abstract base ring for Schubert-family polynomials.

Concrete subclasses supply the multiplication kernels (``double_mul``/``single_mul``,
``mult_poly_single``/``mult_poly_double``), the basis-element expansion
(``cached_schubpoly``), printing (``printing_term``), coercion, and basis changes.
Two rings compare equal iff they have the same type and generating sets.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(genset, coeff_genset, domain=None)
```

**Arguments**:

- `genset` - Primary generating set (the ``x`` variables).
- `coeff_genset` - Coefficient generating set (``y``/``z``), or a set with ``label=None`` for single rings.
- `domain` - Optional coefficient domain passed to `BaseRing`.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via `_mul_schub_dicts`, which dispatches to the appropriate
`schubmult.mult` kernel based on both rings' generating sets; scalars go through `BaseRing.mul`.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.new"></a>

#### new

```python
def new(x)
```

Hook: build an element from ``x`` (permutation, code, or expression).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Hook: the sympy symbol displayed for basis element ``k``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Hook: coproduct of basis element ``k`` in the tensor square ring.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(elem)
```

Hook: whether ``elem`` should be multiplied via the elementary-symmetric fast path.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Hook: elementary-symmetric fast-path multiplication.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.elem_sym"></a>

#### elem\_sym

```python
@property
def elem_sym()
```

Hook: the elementary symmetric polynomial function used by this ring.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

Hook: symbolic (unevaluated) elementary symmetric function for ``in_SEM_basis``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Hook: substitution dict turning the symbolic elementary symmetric symbols back into polynomials.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce ``element`` into the coefficient domain, refusing anything containing a ``genset`` variable.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.genset"></a>

#### genset

```python
@property
def genset()
```

Primary generating set (the ``x`` variables).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.coeff_genset"></a>

#### coeff\_genset

```python
@property
def coeff_genset()
```

Coefficient generating set (``y``/``z`` variables); ``label`` is ``None`` for single rings.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Hook: re-express ``elem`` in the quantum Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Hook: re-express ``elem`` in the classical Schubert basis.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.quantum_schubpoly"></a>

#### quantum\_schubpoly

```python
def quantum_schubpoly(perm)
```

Hook: the quantum Schubert polynomial for ``perm``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.cached_product"></a>

#### cached\_product

```python
def cached_product(u, v, basis2)
```

Hook: cached structure constants of ``S_u * S_v`` with ``v`` in ``basis2``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
def cached_positive_product(u, v, basis2)
```

Hook: like ``cached_product`` but with manifestly positive coefficients.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Hook: multiply ``elem`` by a symbolic expression ``x``.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

Hook: the double-Schubert multiplication kernel (e.g. ``schubmult_double``).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

Hook: the single-Schubert multiplication kernel (e.g. ``schubmult_py``).

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

Hook: the single-variant multiply-by-polynomial kernel.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

Hook: the double-variant multiply-by-polynomial kernel.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.quantum_elem_func"></a>

#### quantum\_elem\_func

```python
@property
def quantum_elem_func()
```

Hook: the quantum elementary symmetric function.

<a id="schubmult.rings.schubert.base_schubert_ring.BaseSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
def cached_schubpoly(k)
```

Hook: the (cached) explicit polynomial for basis element ``k``.

<a id="schubmult.rings.schubert.beta_coxeter"></a>

# schubmult.rings.schubert.beta\_coxeter

Scaffold for a beta-deformed (Grothendieck / 0-Hecke) Coxeter operator ring.

Intended to model the beta-isobaric divided differences ``pi_i = partial_i + beta (x_i partial_i - 1)``
and their relations with the simple reflections (see the commented-out relations above
`BetaCoxeterRing`). At present the implementation is an unmodified copy of
`schubmult.rings.schubert.nil_hecke` -- `BetaCoxeterRing`/`BetaCoxeterElement` behave
identically to `NilHeckeRing`/`NilHeckeElement`, and the deformation has not been wired in.
Prefer `nil_hecke` for actual use; this module is kept as a starting point for that work.

<a id="schubmult.rings.schubert.beta_coxeter.BetaCoxeterElement"></a>

## BetaCoxeterElement Objects

```python
class BetaCoxeterElement(DomainElement, DefaultPrinting, dict)
```

An element of a `BetaCoxeterRing`; currently identical in behavior to `NilHeckeElement`.

<a id="schubmult.rings.schubert.beta_coxeter.BetaCoxeterRing"></a>

## BetaCoxeterRing Objects

```python
class BetaCoxeterRing(Ring, CompositeDomain)
```

Scaffold ring; currently identical in behavior to `NilHeckeRing` (see module docstring).

<a id="schubmult.rings.schubert.chevalley"></a>

# schubmult.rings.schubert.chevalley

Lenart--Postnikov :math:`K_T`-Chevalley formula in type :math:`A`.

Implements Theorem 6.1 / Proposition 14.5 of Lenart--Postnikov, *Affine Weyl
groups in K-theory and representation theory* (arXiv:math/0309207)::

    e^lambda . [O_u] = sum_{w, mu} c^{lambda, mu}_{u, w} x^mu [O_w]

    c^{lambda, mu}_{u, w} = sum_J (-1)^{n(J)}

summed over subsets ``J = {j_1 < ... < j_s}`` of a fixed lambda-chain such that

(a) ``u > u r_{j_1} > ... > u r_{j_1} ... r_{j_s} = w`` is a saturated
    decreasing chain in Bruhat order, and
(b) ``-mu = u r_{j_1} ... r_{j_s} (-lambda)``,

with ``n(J)`` the number of negative roots among ``beta_{j_1}, ..., beta_{j_s}``.

Weights are integer vectors in the ``epsilon`` basis; a root ``(a, b)`` denotes
``eps_a - eps_b`` and is positive iff ``a < b``.

<a id="schubmult.rings.schubert.chevalley.fundamental_weight"></a>

#### fundamental\_weight

```python
def fundamental_weight(k, n)
```

``omega_k = eps_1 + ... + eps_k`` as a length-``n`` vector.

<a id="schubmult.rings.schubert.chevalley.lambda_chain"></a>

#### lambda\_chain

```python
@cache
def lambda_chain(weight, n)
```

Reduced ``lambda``-chain for ``weight`` in ``A_{n-1}``, via Prop. 6.7.

Returns a tuple of ``(beta, root, k)``: ``root`` is the positive root
``alpha`` of the affine reflection ``r_j = s_{alpha, k}``, and ``beta`` is
the signed root ``b(r_j)`` whose sign determines ``(-1)^{n(J)}``.

<a id="schubmult.rings.schubert.chevalley.kt_chevalley_coefficients"></a>

#### kt\_chevalley\_coefficients

```python
def kt_chevalley_coefficients(u, weight, n=None)
```

Chevalley coefficients for ``e^weight . [O_u]``.

Returns ``{(w, mu): coefficient}`` with ``mu`` an integer weight vector.

<a id="schubmult.rings.schubert.double_grothendieck_ring"></a>

# schubmult.rings.schubert.double\_grothendieck\_ring

Double (equivariant K-theoretic) Grothendieck polynomial ring: the ``DGx`` interface.

`DoubleGrothendieckRing` represents ``G_w(x; y)`` with deformation parameter
``beta``. Products go through `schubmult.mult.groth_double.grothmult_double`
(the K-theoretic Monk/Pieri machinery) where available, with a fallback that
expands into the underlying `DoubleSchubertRing`. The ring also exposes the
localization/vanishing data used to convert between the Schubert and
Grothendieck bases (``permuted_subs_dict``, ``product_of_roots``,
``exp_root``, ``chevalley``).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement"></a>

## DoubleGrothendieckElement Objects

```python
class DoubleGrothendieckElement(BaseSchubertElement)
```

Element of a DoubleGrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Expand to an explicit polynomial: ``sum coeff * G_w(x; y)``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement.perm_subs"></a>

#### perm\_subs

```python
def perm_subs(perm)
```

Localize at the torus fixed point ``perm``: substitute ``x_i -> (-) y_{perm(i)}`` (formal inverse).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing"></a>

## DoubleGrothendieckRing Objects

```python
class DoubleGrothendieckRing(BaseSchubertRing)
```

Ring of double (K-theoretic) Grothendieck polynomials G_w(x, y).

Elements are stored in the G-basis as ``{Permutation: coeff}``. There is
no direct structure-constant formula for the product implemented here:
instead both factors are expanded into the underlying ``DoubleSchubertRing``
(via ``grothendieck_poly_with_ring``), multiplied there, and the product is
converted back to the G-basis with ``to_groth_with_ring``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.perm_subs"></a>

#### perm\_subs

```python
def perm_subs(elem, perm)
```

Localize ``elem`` at ``perm``: expand into double Schubert polynomials and substitute
``x_i -> -y_{perm(i)} / (1 + beta y_{perm(i)})``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.groth_double.grothmult_double`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.groth.grothmult_py`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.beta"></a>

#### beta

```python
@property
def beta()
```

The deformation parameter.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.single_variable"></a>

#### single\_variable

```python
@property
def single_variable()
```

`schubmult.mult.groth_double.single_variable_groth` (K-theoretic Monk rule for ``x_k``).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.vanish_subs_dict"></a>

#### vanish\_subs\_dict

```python
@cached_property
def vanish_subs_dict()
```

Localization at the identity: ``x_i -> -y_i / (1 + beta y_i)`` for the first 50 variables.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.permuted_subs_dict"></a>

#### permuted\_subs\_dict

```python
def permuted_subs_dict(perm, length=None)
```

Localization at ``perm``: ``x_i -> -y_{perm(i)} / (1 + beta y_{perm(i)})`` for ``i <= length``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.product_of_roots"></a>

#### product\_of\_roots

```python
@cache
def product_of_roots(perm)
```

``prod_{(a,b) in Inv(perm^-1)} (y_a (-) y_b)``: the localization of ``G_perm`` at itself (Euler class).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_root"></a>

#### exp\_root

```python
@cache
def exp_root(a, b)
```

Polynomial form of the root ``eps_a - eps_b``, i.e. ``(x^root - 1)/beta``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_weight"></a>

#### exp\_weight

```python
@cache
def exp_weight(index)
```

Polynomial form of the fundamental weight ``eps_index``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_character"></a>

#### exp\_character

```python
def exp_character(weight)
```

Polynomial form of the multiplicative character ``x^weight``.

Determined by ``exp_root``: ``x^{eps_i} = 1 + beta * y_i``, so that
``(x^{eps_a - eps_b} - 1)/beta`` reproduces ``exp_root(a, b)``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.chevalley"></a>

#### chevalley

```python
def chevalley(weight, perm, n=None)
```

Lenart--Postnikov K_T-Chevalley formula for ``e^weight * G_perm``.

``weight`` is an integer vector in the ``epsilon`` basis.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.div_by_product_of_roots"></a>

#### div\_by\_product\_of\_roots

```python
def div_by_product_of_roots(expr, perm)
```

Divide ``expr`` by ``product_of_roots(perm)`` one linear factor at a time, leaving any
non-exact factors in the denominator (so the result may be a rational function).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.groth_double.mult_poly_groth_double`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.groth.mult_poly_groth`.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.schub_as_groth"></a>

#### schub\_as\_groth

```python
@cache
def schub_as_groth(perm)
```

The double Schubert polynomial ``S_perm`` expanded in the Grothendieck basis (cached).

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.from_double_schubert_elem"></a>

#### from\_double\_schubert\_elem

```python
def from_double_schubert_elem(elem)
```

Convert a `DoubleSchubertElement` into this ring's Grothendieck basis.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression: single ``x`` variables via the K-Monk rule, ``Add``/``Mul``/``Pow``
recursively, anything else as a coefficient.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Ring product via ``_best_effort_grothmult_double``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial into the Grothendieck basis by multiplying the identity by it.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit ``G_k(x; y)``, as a sum of `WCGraph` monomials weighted by ``beta^excess``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="")
```

The ``DoubleGrothendieckPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an element of this ring, or a polynomial expression.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.from_dict"></a>

#### from\_dict

```python
def from_dict(dct)
```

Build an element from ``{Permutation: coeff}``, dropping exact zeros.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DGx"></a>

#### DGx

```python
def DGx(x, genset=GeneratingSet("y"))
```

Construct a double Grothendieck element in ``x`` with coefficient alphabet ``genset`` (a
`GeneratingSet`, a label string, or ``"0"`` for the zero alphabet).

<a id="schubmult.rings.schubert.double_schubert_ring"></a>

# schubmult.rings.schubert.double\_schubert\_ring

Double Schubert polynomial ring: the ``DSx`` interface.

`DoubleSchubertRing` represents ``Z[y][x]`` in the basis of double Schubert
polynomials ``S_w(x; y)``, dispatching products to `schubmult.mult.double`.
It is also the workhorse behind the single ring (`schubert_ring.SingleSchubertRing`
is a `DoubleSchubertRing` with an all-zero coefficient alphabet). Beyond ring
arithmetic, `DoubleSchubertElement` supports divided differences, isobaric
divided differences, variable substitution/evaluation, coproducts, and
expansion into elementary-symmetric ("CEM"/"SEM") bases.

Variants: `ElemDoubleSchubertRing` keeps coefficients as unevaluated factorial
elementary symmetric functions; `DoubleSchubertRingDown` uses the descent-side
("down") kernels.

<a id="schubmult.rings.schubert.double_schubert_ring.is_fact_elem_sym"></a>

#### is\_fact\_elem\_sym

```python
def is_fact_elem_sym(obj)
```

Whether ``obj`` is an (unevaluated) factorial elementary symmetric function.

<a id="schubmult.rings.schubert.double_schubert_ring.is_fact_complete_sym"></a>

#### is\_fact\_complete\_sym

```python
def is_fact_complete_sym(obj)
```

Whether ``obj`` is an (unevaluated) factorial complete homogeneous symmetric function.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement"></a>

## DoubleSchubertElement Objects

```python
class DoubleSchubertElement(BaseSchubertElement)
```

An element of a `DoubleSchubertRing`: ``{Permutation: coefficient}`` in the
double Schubert basis ``S_w(x; y)``, with sympy coefficients in ``y``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.to_genset_dict"></a>

#### to\_genset\_dict

```python
def to_genset_dict(trim=False)
```

Expand to a polynomial and return ``{exponent_tuple: coeff}`` over the ``x`` variables;
``trim=True`` merges keys that differ only by trailing zeros.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.divdiff"></a>

#### divdiff

```python
def divdiff(i)
```

Divided difference ``partial_i``: ``S_w -> S_{w s_i}`` when ``i`` is a descent of ``w``, else 0.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.simpleref"></a>

#### simpleref

```python
def simpleref(i)
```

Action of the simple reflection ``s_i`` on the ``x`` variables: ``f + (x_{i+1} - x_i) partial_i f``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.coeff_isobaric"></a>

#### coeff\_isobaric

```python
def coeff_isobaric(i, beta)
```

Isobaric divided difference acting on the ``y`` (coefficient) alphabet, transported through
the basis via the antipode-style inversion ``S_w -> (-1)^{l(w)} S_{w^{-1}}``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.isobaric"></a>

#### isobaric

```python
def isobaric(i, beta)
```

Beta-deformed isobaric divided difference ``pi_i = partial_i + beta (x_i partial_i - 1)``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply ``partial_w`` for ``w = perm``, peeling simple reflections from the last descent.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.isobaric_perm"></a>

#### isobaric\_perm

```python
def isobaric_perm(perm, beta)
```

Apply the beta-isobaric ``pi_w`` for ``w = perm``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.isobaric_plus_beta"></a>

#### isobaric\_plus\_beta

```python
def isobaric_plus_beta(i, beta)
```

The variant ``partial_i + beta x_i partial_i`` (isobaric without the ``-beta`` identity term).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.act"></a>

#### act

```python
def act(perm)
```

Permute the ``x`` variables by ``perm``, as a composition of ``simpleref``s.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.max_index"></a>

#### max\_index

```python
def max_index()
```

The largest ``x`` index (1-indexed) any basis permutation actually depends on.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.eval"></a>

#### eval

```python
def eval(x)
```

Substitute ``{generator: value}`` pairs one at a time (via ``pull_out_gen``); returns a
scalar if the result collapses to the identity basis element.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.subs"></a>

#### subs

```python
def subs(old, new)
```

Substitute ``old -> new`` where ``old`` is an ``x`` variable (moved to the last position and
pulled out via ``pull_out_var``), a ``y`` variable (transported through the basis), or a plain
coefficient symbol.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.free_symbols"></a>

#### free\_symbols

```python
@property
def free_symbols()
```

Coefficient symbols plus the ``x``/``y`` variables the basis permutations actually depend on.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.pull_out_gen"></a>

#### pull\_out\_gen

```python
def pull_out_gen(gen)
```

Factor out all dependence on one generator ``gen`` (an ``x`` or ``y`` variable), returning an
element over a `MaskedGeneratingSet` ring with ``gen`` removed and explicit ``(gen - y_j)``
(or factorial-elementary-symmetric) prefactors.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.in_CEM_basis"></a>

#### in\_CEM\_basis

```python
def in_CEM_basis()
```

Expand in the complete-elementary-monomial (CEM) basis using the ring's symbolic elementary function.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.cem_rep"></a>

#### cem\_rep

```python
def cem_rep(elem_func, mumu=None)
```

CEM expansion with a custom ``elem_func``; ``mumu`` selects a dominant permutation to expand
against (defaults to the classical route).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.coproduct"></a>

#### coproduct

```python
def coproduct(*indices,
              alt_coeff_genset=None,
              on_coeff_gens=False,
              gname1=None,
              gname2=None)
```

Coproduct splitting the ``x`` variables (or ``y`` if ``on_coeff_gens``) at the given 1-indexed
``indices``: returns an element of the `TensorRing` of two `DoubleSchubertRing`s over the
complementary `MaskedGeneratingSet`s, labeled ``gname1``/``gname2``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.max_gens"></a>

#### max\_gens

```python
@cached_property
def max_gens()
```

Largest 0-indexed descent over all basis permutations.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.positive_elem_sym_rep"></a>

#### positive\_elem\_sym\_rep

```python
def positive_elem_sym_rep()
```

Manifestly positive expansion in factorial elementary symmetric functions (forward ``pull_out_var``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.positive_elem_sym_rep_backward"></a>

#### positive\_elem\_sym\_rep\_backward

```python
def positive_elem_sym_rep_backward()
```

Like ``positive_elem_sym_rep`` but peeling from the last descent backward.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement.antipode"></a>

#### antipode

```python
def antipode()
```

The antipode: swap the two alphabets and invert each basis permutation (see `DoubleSchubertRing.antipode`).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing"></a>

## DoubleSchubertRing Objects

```python
class DoubleSchubertRing(BaseSchubertRing)
```

The ring of double Schubert polynomials ``S_w(x; y)`` over ``genset`` (``x``) and
``coeff_genset`` (``y``). Call the ring with a permutation, Lehmer code, or polynomial
expression to construct an element; the module-level ``DSx`` is the standard instance.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.coeff_ring"></a>

#### coeff\_ring

```python
@cached_property
def coeff_ring()
```

The single Schubert ring over the coefficient alphabet ``y`` (used by ``coeff_isobaric``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.antipode_ring"></a>

#### antipode\_ring

```python
@cached_property
def antipode_ring()
```

The same ring with the two alphabets swapped.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.antipode"></a>

#### antipode

```python
def antipode(elem)
```

Map ``sum c_w S_w(x; y)`` to ``sum c_w S_{w^{-1}}(y; x)`` in the swapped-alphabet ring.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Right-multiply by a scalar (coefficient-domain element) or, failing that, by an expression.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.positive_elem_sym_rep"></a>

#### positive\_elem\_sym\_rep

```python
def positive_elem_sym_rep(perm, index=1)
```

Manifestly positive expansion of ``S_perm`` in factorial elementary symmetric functions, peeling
the first variable of ``~perm`` at each step (``pull_out_var(1, ...)``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.positive_elem_sym_rep_backward"></a>

#### positive\_elem\_sym\_rep\_backward

```python
def positive_elem_sym_rep_backward(perm)
```

Like ``positive_elem_sym_rep`` but peeling from the last descent of ``~perm`` backward.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="")
```

The ``DSchubPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
@property
def elem_sym()
```

`FactorialElemSym`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(other)
```

Whether ``other`` is a factorial elementary symmetric function (eligible for ``elem_mul``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Multiply by a factorial elementary symmetric function in ``x`` variables via the positional
Pieri rule (``elem_sym_positional_perms``), expanding the leftover factor with ``expand_func``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

`FactorialElemSym` (kept unevaluated for symbolic expansions).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.schubert_schur_elem_func"></a>

#### schubert\_schur\_elem\_func

```python
def schubert_schur_elem_func(numvars)
```

Elementary-symmetric substitute for the Schubert-tensor-Schur expansion: ``e_p(x_1..x_k)`` maps
to a Schubert basis element on the left factor when ``k >= numvars`` and on the right otherwise.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_schubert_schur_basis"></a>

#### in\_schubert\_schur\_basis

```python
def in_schubert_schur_basis(perm, numvars)
```

Expand ``S_perm`` in the Schubert-tensor-Schur basis, treating the last ``numvars`` variables
as the symmetric (Schur) part.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_descending_schur_basis"></a>

#### in\_descending\_schur\_basis

```python
def in_descending_schur_basis(perm, numvars)
```

Iterate ``in_schubert_schur_basis`` down through ``numvars, numvars-1, ..., 1``, producing a
nested tensor of Schur-like factors.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Substitution dict ``{e_p_k: elem_sym_poly(p, k, x)}`` for all ``1 <= p <= k <= kk``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.flip"></a>

#### flip

```python
@staticmethod
def flip(elem)
```

Re-express a factorial elementary symmetric function with its two alphabets swapped, via the
corresponding Grassmannian Schubert polynomial's CEM expansion.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Expand each basis element via ``quantum_schubpoly`` (a quantum double Schubert element).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Identity (this ring is already classical).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.quantum_schubpoly"></a>

#### quantum\_schubpoly

```python
@cache
def quantum_schubpoly(perm)
```

The classical ``S_perm`` expressed in the quantum double Schubert basis (via ``quantum_elem_func``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants of ``S_u(x; y) * S_v(x; z)`` (``z`` = ``basis2.coeff_genset``), via ``schubmult_double``.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Like ``cached_product`` but with manifestly positive coefficients (generic alphabets, then substituted).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.double.schubmult_double`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.single.schubmult_py`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.single.mult_poly_py`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.double.mult_poly_double`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.quantum_elem_func"></a>

#### quantum\_elem\_func

```python
@property
def quantum_elem_func()
```

Elementary symmetric function valued in the quantum double Schubert ring, computed by a
divide-and-conquer recursion on the variable set (used by ``quantum_schubpoly``).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.monomial_schub"></a>

#### monomial\_schub

```python
def monomial_schub(monom)
```

The monomial ``x^monom`` expressed in the Schubert basis (trailing zeros in ``monom`` ignored).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit polynomial ``S_k(x; y)`` (cached).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.complete_mul"></a>

#### complete\_mul

```python
def complete_mul(elem, x)
```

Multiply by a factorial complete homogeneous symmetric function in ``x`` variables via
``complete_sym_positional_perms`` (the dual Pieri rule).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.handle_sympoly"></a>

#### handle\_sympoly

```python
def handle_sympoly(other)
```

How a symmetric-function coefficient is stored: evaluated to a polynomial here.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.single_variable"></a>

#### single\_variable

```python
def single_variable(elem, varnum)
```

Multiply by the single variable ``x_varnum`` (equivariant Monk rule).

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial expression in ``x``/``y`` into the Schubert basis.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply ``elem`` by an arbitrary expression ``x``: single variables use the Monk rule,
(factorial) elementary/complete symmetric functions use their Pieri rules (splitting out
variables from the wrong alphabet as needed), and ``Add``/``Mul``/``Pow`` recurse; anything
else is treated as a coefficient.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an existing element of this ring, or an expression.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown"></a>

## DoubleSchubertRingDown Objects

```python
class DoubleSchubertRingDown(DoubleSchubertRing)
```

`DoubleSchubertRing` using the descent-side ("down") multiplication kernels
(``schubmult_double_down``/``schubmult_py_down``); basis symbols print with an ``op`` prefix.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.double.schubmult_double_down`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.single.schubmult_py_down`.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Down-kernel structure constants over generic alphabets, substituted back to the ring's alphabets.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Positive variant of ``cached_product`` for the down kernel.

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertRingDown.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="op")
```

The ``DSchubPoly`` display symbol, prefixed with ``op`` by default.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing"></a>

## ElemDoubleSchubertRing Objects

```python
class ElemDoubleSchubertRing(DoubleSchubertRing)
```

`DoubleSchubertRing` whose coefficients are kept as unevaluated `FactorialElemSym`
functions instead of being expanded to polynomials; products use the ``*_from_elems`` kernels.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.replacematch"></a>

#### replacematch

```python
@property
def replacematch()
```

A ``(a, b) -> expression`` rewriter turning differences ``a - b`` into `FactorialElemSym(1, 1, ...)`
forms, respecting which alphabet each symbol belongs to.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.elem_func"></a>

#### elem\_func

```python
@property
def elem_func()
```

`FactorialElemSym`.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.handle_sympoly"></a>

#### handle\_sympoly

```python
def handle_sympoly(other)
```

Keep symmetric-function coefficients unevaluated.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Positional Pieri rule for a factorial elementary symmetric function, keeping the leftover
factor as an unevaluated coefficient.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.complete_mul"></a>

#### complete\_mul

```python
def complete_mul(elem, x)
```

Dual Pieri rule for a factorial complete symmetric function, keeping the leftover factor unevaluated.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants via ``schubmult_double_from_elems`` with `FactorialElemSym` coefficients.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Structure constants via the positive ``schubmult_double_alt_from_elems`` route.

<a id="schubmult.rings.schubert.double_schubert_ring.ElemDoubleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an element of this ring, or an expression.

<a id="schubmult.rings.schubert.double_schubert_ring.DSx"></a>

#### DSx

```python
def DSx(x, genset=GeneratingSet("y"), elem_sym=False, down=False)
```

Construct a double Schubert polynomial element in ``x`` with coefficient alphabet ``genset``.

``DSx([3, 1, 2])`` is ``S_{312}(x; y)``. Pass ``genset="z"`` (or a `GeneratingSet`) for a
different coefficient alphabet; ``elem_sym=True`` uses `ElemDoubleSchubertRing`, ``down=True``
uses `DoubleSchubertRingDown`.

<a id="schubmult.rings.schubert.grothendieck_ring"></a>

# schubmult.rings.schubert.grothendieck\_ring

Grothendieck polynomial ring (non-equivariant): the ``Gx`` interface.

`GrothendieckRing` is the beta-deformation of `SingleSchubertRing`; basis
elements ``G_w`` are the K-theoretic Schubert classes, with ``beta = 0``
recovering ordinary Schubert polynomials. There is no coefficient alphabet yet
(see `double_grothendieck_ring` for the equivariant version).

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement"></a>

## GrothendieckElement Objects

```python
class GrothendieckElement(BaseSchubertElement)
```

Element of a GrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement.as_polynomial"></a>

#### as\_polynomial

```python
def as_polynomial()
```

Expand to an explicit polynomial: ``sum coeff * G_w(x)``.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement.mult_poly"></a>

#### mult\_poly

```python
def mult_poly(poly)
```

Multiply by an arbitrary polynomial in ``x`` via the Grothendieck Chevalley rule (`mult_poly_groth`).

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing"></a>

## GrothendieckRing Objects

```python
class GrothendieckRing(BaseSchubertRing)
```

Ring of Grothendieck polynomials.

A deformation of the Schubert polynomial ring with parameter beta.
Basis elements G_w satisfy G_u * G_v = sum_w c^w_{u,v}(beta) G_w
where c^w_{u,v}(beta) are polynomials in beta with integer coefficients.
When beta=0, recovers ordinary Schubert polynomials.

Parameters
----------
genset : GeneratingSet
    The generating set (variable alphabet).
beta : sympy/symengine symbol, optional
    The deformation parameter. Defaults to Symbol("β").

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.beta"></a>

#### beta

```python
@property
def beta()
```

The deformation parameter.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`mult_poly_groth` with this ring's ``beta`` bound.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial to the Grothendieck basis: expand in Schubert polynomials first, then
change basis Schubert -> Grothendieck.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, expr)
```

Multiply by an expression by first converting it into the Grothendieck basis.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants ``c^w_{u,v}(beta)`` via ``groth_mul_full_with_ring``; only same-ring products supported.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product``.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit Grothendieck polynomial ``G_k(x)`` (cached).

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k, prefix="")
```

The ``GrothendieckPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, an element of this ring, or a polynomial expression.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing.from_dict"></a>

#### from\_dict

```python
def from_dict(dct)
```

Build an element from ``{Permutation: coeff}``, dropping terms whose coefficient expands to zero.

<a id="schubmult.rings.schubert.nil_hecke"></a>

# schubmult.rings.schubert.nil\_hecke

The nilHecke ring of divided-difference operators acting on Schubert polynomials.

`NilHeckeRing` elements are ``{Permutation: coefficient}`` combinations of the
divided-difference operators ``partial_w`` (printed ``df(w)``), with polynomial
coefficients in the ``x`` variables multiplied on the left. ``partial_w`` acts on
a `DoubleSchubertElement` via `NilHeckeElement.apply`, sending ``S_v -> S_{v w^{-1}}``
when length-additive. Products use the descent-side kernel ``schubmult_double_down``
to commute polynomial coefficients past operators. The module-level ``df`` is the
standard instance in ``x``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement"></a>

## NilHeckeElement Objects

```python
class NilHeckeElement(DomainElement, DefaultPrinting, dict)
```

An element of a `NilHeckeRing`: ``{Permutation: coeff}`` combination of divided-difference operators.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.apply"></a>

#### apply

```python
def apply(other)
```

Act on a `DoubleSchubertElement`: each ``partial_w`` sends ``S_v -> S_{v w^{-1}}`` when
``l(v w^{-1}) = l(v) - l(w)``, else kills it; coefficients multiply the result.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_terms"></a>

#### as\_terms

```python
def as_terms()
```

Terms ``coeff * df(w)`` in dict order (sympy printing hook).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by permutation length then lexicographically (sympy printing hook).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

``{df(w): coeff}`` mapping display symbols to coefficients.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand each coefficient, keeping the operator basis.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeElement.as_expr"></a>

#### as\_expr

```python
def as_expr()
```

Sum of the ``as_terms()`` as a sympy ``Add``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing"></a>

## NilHeckeRing Objects

```python
class NilHeckeRing(Ring, CompositeDomain)
```

The nilHecke ring in the alphabet ``genset``; see the module docstring. ``df`` is the standard instance.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.to_sympy"></a>

#### to\_sympy

```python
def to_sympy(elem)
```

Convert an element to a sympy expression (``as_expr``).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.isobaric"></a>

#### isobaric

```python
def isobaric(perm, groth=False, *, groth_beta=None)
```

The isobaric divided difference ``pi_perm`` as a nilHecke element: ``pi_i = partial_i x_{i+1}``
(or the Grothendieck version ``partial_i (1 + beta x_{i+1})`` with ``groth=True``), composed
along a reduced word of ``perm``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.g_isobaric"></a>

#### g\_isobaric

```python
def g_isobaric(perm)
```

``isobaric(perm, groth=True)``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.fgp_operator"></a>

#### fgp\_operator

```python
def fgp_operator(k, length, q_var=GeneratingSet("q"))
```

The Fomin-Gelfand-Postnikov quantization of ``x_k`` as a nilHecke element:
``x_k - sum_{i<k} q_i...q_{k-1} partial_{(i k)} + sum_{i>k} q_k...q_{i-1} partial_{(k i)}``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.subs_fgp"></a>

#### subs\_fgp

```python
def subs_fgp(poly, length)
```

Substitute every ``x_k`` in ``poly`` by its ``fgp_operator`` (quantize a polynomial).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.mul_scalar"></a>

#### mul\_scalar

```python
def mul_scalar(elem, other)
```

Multiply on the right by a polynomial/Schubert element, commuting it past the operators via
``schubmult_double_down`` (Leibniz rule for divided differences).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.mul_perm"></a>

#### mul\_perm

```python
def mul_perm(elem, perm)
```

Right-multiply every operator ``partial_k`` by ``partial_perm``, keeping only length-additive products.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.rmul"></a>

#### rmul

```python
def rmul(elem, other)
```

Left-multiply by a scalar/polynomial (coefficients sit on the left, so this is plain scaling).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.mul"></a>

#### mul

```python
def mul(elem, other)
```

Ring product: scalars scale, nilHecke elements combine via ``mul_scalar`` then ``mul_perm``, else ``mul_scalar``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list (the operator ``partial_w``) or a polynomial (a scalar).

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The display symbol ``df(w)`` / ``∂(w)`` / ``\partial^w`` for the operator indexed by ``k``.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce ``element`` into the coefficient domain, refusing ring elements and anything containing an ``x`` variable.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.genset"></a>

#### genset

```python
@property
def genset()
```

The ``x`` alphabet.

<a id="schubmult.rings.schubert.nil_hecke.NilHeckeRing.from_expr"></a>

#### from\_expr

```python
def from_expr(x)
```

Build the scalar element ``x * partial_id``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring"></a>

# schubmult.rings.schubert.parabolic\_quantum\_double\_schubert\_ring

Parabolic quantum double Schubert polynomial ring: the ``QPDSx`` interface.

`ParabolicQuantumDoubleSchubertRing` models the quantum cohomology of a partial
flag variety with block sizes ``index_comp``. Basis permutations must be
parabolic (increasing within each block). Products are computed in the full
flag quantum ring and projected down via the Peterson-Woodward comparison
(`schubmult.mult.quantum_double.apply_peterson_woodward`, through
``process_coeff_dict``). The parabolic quantum elementary symmetric functions
acquire a ``q``-correction at each block boundary.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertElement"></a>

## ParabolicQuantumDoubleSchubertElement Objects

```python
class ParabolicQuantumDoubleSchubertElement(BaseSchubertElement)
```

An element of a `ParabolicQuantumDoubleSchubertRing`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertElement.index_comp"></a>

#### index\_comp

```python
@property
def index_comp()
```

The ring's block-size composition.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertElement.kill_ideal"></a>

#### kill\_ideal

```python
def kill_ideal()
```

Drop basis permutations longer than ``sum(index_comp)`` (those lie in the ideal cut out by the parabolic).

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing"></a>

## ParabolicQuantumDoubleSchubertRing Objects

```python
class ParabolicQuantumDoubleSchubertRing(BaseSchubertRing)
```

Quantum double Schubert polynomials for the partial flag variety with block sizes ``index_comp``.
Construct via ``QPDSx(*index_comp)([perm])`` or ``make_parabolic_quantum_basis``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(genset, coeff_genset, index_comp)
```

**Arguments**:

- `genset` - Primary ``x`` alphabet.
- `coeff_genset` - Coefficient ``y`` alphabet.
- `index_comp` - Composition of block sizes; ``sum(index_comp)`` is the ambient ``n``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

Symbolic elementary symmetric function ``e_p_k`` combined with complete symmetric corrections in ``-y``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Substitution dict ``{e_p_k: elem_sym(p, k, x, 0)}`` for all ``1 <= p <= k <= kk``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.parabolic_index"></a>

#### parabolic\_index

```python
@property
def parabolic_index()
```

1-indexed positions of the simple reflections inside the parabolic subgroup (within-block positions).

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.quantum_basis"></a>

#### quantum\_basis

```python
@property
def quantum_basis()
```

The full-flag `QuantumDoubleSchubertRing` over the same alphabets.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.classical_basis"></a>

#### classical\_basis

```python
@property
def classical_basis()
```

The classical `DoubleSchubertRing` over the same alphabets.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
def elem_sym(p, k, varl1, varl2)
```

Parabolic quantum double elementary symmetric polynomial ``E_p(x_1..x_k; y)``: classical below
the first block boundary, with a ``q_j``-correction at each block boundary ``N_j``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.index_comp"></a>

#### index\_comp

```python
@property
def index_comp()
```

The block-size composition.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.process_coeff_dict"></a>

#### process\_coeff\_dict

```python
def process_coeff_dict(coeff_dict)
```

Project a full-flag quantum coefficient dict onto this parabolic ring via Peterson-Woodward,
extending the parabolic index if any permutation exceeds the ambient ``n``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Full-flag quantum double product (generic alphabets, substituted) then projected via ``process_coeff_dict``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Expand into the full-flag quantum double Schubert basis via ``quantum_elem_func``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Expand into the classical double Schubert basis via ``quantum_as_classical_schubpoly``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.classical_in_basis"></a>

#### classical\_in\_basis

```python
@cache
def classical_in_basis(k)
```

Express the classical ``S_k`` in this parabolic quantum basis, by iteratively subtracting
off lower-order corrections until the polynomials agree.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.classical_elem_func"></a>

#### classical\_elem\_func

```python
@property
def classical_elem_func()
```

Parabolic quantum elementary symmetric function valued in the classical `DoubleSchubertRing`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.quantum_elem_func"></a>

#### quantum\_elem\_func

```python
@property
def quantum_elem_func()
```

Parabolic quantum elementary symmetric function valued in the full-flag `QuantumDoubleSchubertRing`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The ``PQDSchubPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.quantum_as_classical_schubpoly"></a>

#### quantum\_as\_classical\_schubpoly

```python
@cache
def quantum_as_classical_schubpoly(perm)
```

``S^{q,P}_perm`` expanded in the classical double Schubert basis, against the appropriate longest element.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit parabolic quantum double Schubert polynomial for ``k``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Positive variant of ``cached_product`` via ``schubmult_q_generic_partial_posify``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

``schubmult_q_double_fast`` followed by ``process_coeff_dict``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

``schubmult_q_fast`` followed by ``process_coeff_dict``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.quantum.mult_poly_q`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.quantum_double.mult_poly_q_double`.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial to this basis by peeling off leading monomials; raises ``ValueError`` if
``expr`` lacks the within-block symmetry the parabolic ring requires.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression by first converting it into this basis.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.ParabolicQuantumDoubleSchubertRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

Build an element from a parabolic permutation/Lehmer list or an expression; raises ``ValueError``
if the permutation is not parabolic for this ring's blocks.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.make_parabolic_quantum_basis"></a>

#### make\_parabolic\_quantum\_basis

```python
def make_parabolic_quantum_basis(index_comp, coeff_genset)
```

The `ParabolicQuantumDoubleSchubertRing` in ``x`` for block sizes ``index_comp`` and coefficient alphabet ``coeff_genset``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.QPDSx_index"></a>

#### QPDSx\_index

```python
def QPDSx_index(*args)
```

Return a constructor ``f(x, coeff_genset="y")`` building parabolic quantum double elements for block sizes ``args``.

<a id="schubmult.rings.schubert.parabolic_quantum_double_schubert_ring.QPDSx"></a>

#### QPDSx

```python
@cache
def QPDSx(*args)
```

Cached constructor for block sizes ``args``; e.g. ``QPDSx(2, 1)([2, 1, 3])``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring"></a>

# schubmult.rings.schubert.parabolic\_quantum\_schubert\_ring

Parabolic quantum (single) Schubert polynomial ring: the ``QPSx`` interface.

`ParabolicQuantumSingleSchubertRing` is a `ParabolicQuantumDoubleSchubertRing`
with a zero coefficient alphabet, indexed by a composition ``index_comp``
specifying the parabolic subgroup's block sizes.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing"></a>

## ParabolicQuantumSingleSchubertRing Objects

```python
class ParabolicQuantumSingleSchubertRing(ParabolicQuantumDoubleSchubertRing)
```

Quantum Schubert polynomials for the partial flag variety with block sizes ``index_comp``.
Construct via ``QPSx(*index_comp)``; basis permutations must be parabolic for the given blocks.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit parabolic quantum Schubert polynomial for ``k``, expanded against the appropriate
longest element (extended if ``k`` is larger than the ring's default).

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
def elem_sym(p, k, varl1, varl2)
```

Parabolic quantum elementary symmetric polynomial ``E_p(x_1..x_k)``: classical below the first
block boundary, with a ``q``-correction term at each block boundary ``N_j``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.coeff_genset"></a>

#### coeff\_genset

```python
@property
def coeff_genset()
```

Always the zero alphabet.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Full-flag quantum product (``schubmult_q_fast`` or the generic double kernel), then projected
to the parabolic ring via ``process_coeff_dict`` (Peterson-Woodward).

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.ParabolicQuantumSingleSchubertRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

Build an element from a parabolic permutation/Lehmer list or an expression; raises ``ValueError``
if the permutation is not parabolic for this ring's blocks.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.make_single_parabolic_quantum_basis"></a>

#### make\_single\_parabolic\_quantum\_basis

```python
def make_single_parabolic_quantum_basis(index_comp)
```

The `ParabolicQuantumSingleSchubertRing` in ``x`` for block sizes ``index_comp``.

<a id="schubmult.rings.schubert.parabolic_quantum_schubert_ring.QPSx"></a>

#### QPSx

```python
@cache
def QPSx(*args)
```

Cached `ParabolicQuantumSingleSchubertRing` for block sizes ``args``, e.g. ``QPSx(2, 1)([2, 1, 3])``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring"></a>

# schubmult.rings.schubert.quantum\_double\_schubert\_ring

Quantum double Schubert polynomial ring: the ``QDSx`` interface.

`QuantumDoubleSchubertRing` represents quantum double Schubert polynomials
``S^q_w(x; y)`` with quantum parameters ``q_1, q_2, ...`` (the module-level
``q_var``), dispatching products to `schubmult.mult.quantum_double`. Elements
can be converted to/from the classical basis via ``as_classical``/``as_quantum``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.is_fact_elem_sym"></a>

#### is\_fact\_elem\_sym

```python
def is_fact_elem_sym(obj)
```

Whether ``obj`` is an (unevaluated) quantum factorial elementary symmetric function.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertElement"></a>

## QuantumDoubleSchubertElement Objects

```python
class QuantumDoubleSchubertElement(BaseSchubertElement)
```

An element of a `QuantumDoubleSchubertRing`: ``{Permutation: coeff}`` in the basis ``S^q_w(x; y)``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertElement.subs"></a>

#### subs

```python
def subs(old, new)
```

Substitute by round-tripping through the classical basis.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing"></a>

## QuantumDoubleSchubertRing Objects

```python
class QuantumDoubleSchubertRing(BaseSchubertRing)
```

The ring of quantum double Schubert polynomials; ``QDSx`` is the standard constructor.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The ``QDSchubPoly`` display symbol for basis element ``k``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.symbol_elem_func"></a>

#### symbol\_elem\_func

```python
@property
def symbol_elem_func()
```

`QFactorialElemSym`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.elem_sym_subs"></a>

#### elem\_sym\_subs

```python
def elem_sym_subs(kk)
```

Substitution dict ``{e_p_k: elem_sym_poly_q(p, k, x)}`` for all ``1 <= p <= k <= kk``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.elem_sym"></a>

#### elem\_sym

```python
@property
def elem_sym()
```

`QFactorialElemSym`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.is_elem_mul_type"></a>

#### is\_elem\_mul\_type

```python
def is_elem_mul_type(other)
```

Whether ``other`` is a quantum factorial elementary symmetric function.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.elem_mul"></a>

#### elem\_mul

```python
def elem_mul(ring_elem, elem)
```

Multiply by a quantum factorial elementary symmetric function via the quantum positional
Pieri rule (``elem_sym_positional_perms_q``), each term carrying its ``q``-monomial.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants via ``schubmult_q_double_fast``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.in_quantum_basis"></a>

#### in\_quantum\_basis

```python
def in_quantum_basis(elem)
```

Identity (this ring is already quantum).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.in_classical_basis"></a>

#### in\_classical\_basis

```python
def in_classical_basis(elem)
```

Expand each ``S^q_w`` as a classical double Schubert element via ``quantum_as_classical_schubpoly``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.classical_elem_func"></a>

#### classical\_elem\_func

```python
@property
def classical_elem_func()
```

Quantum elementary symmetric function valued in the classical `DoubleSchubertRing`, via the
recursion ``E_p(k) = (x_k - y_{k-p+1}) E_{p-1}(k-1) + E_p(k-1) + q_{k-1} E_{p-2}(k-2)``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.quantum_as_classical_schubpoly"></a>

#### quantum\_as\_classical\_schubpoly

```python
@cache
def quantum_as_classical_schubpoly(perm)
```

``S^q_perm`` expanded in the classical double Schubert basis (cached).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

The explicit quantum double Schubert polynomial ``S^q_k(x; y)`` (cached).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Structure constants with (partially) positive coefficients via ``schubmult_q_generic_partial_posify``.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.double_mul"></a>

#### double\_mul

```python
@property
def double_mul()
```

`schubmult.mult.quantum_double.schubmult_q_double_fast`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.single_mul"></a>

#### single\_mul

```python
@property
def single_mul()
```

`schubmult.mult.quantum.schubmult_q_fast`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.mult_poly_single"></a>

#### mult\_poly\_single

```python
@property
def mult_poly_single()
```

`schubmult.mult.quantum.mult_poly_q`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.positive_elem_sym_rep"></a>

#### positive\_elem\_sym\_rep

```python
def positive_elem_sym_rep(perm, index=1)
```

Manifestly positive expansion of ``S^q_perm`` in quantum factorial elementary symmetric functions (forward).

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.positive_elem_sym_rep_backward"></a>

#### positive\_elem\_sym\_rep\_backward

```python
def positive_elem_sym_rep_backward(perm)
```

Like ``positive_elem_sym_rep`` but peeling from the last descent backward.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.mult_poly_double"></a>

#### mult\_poly\_double

```python
@property
def mult_poly_double()
```

`schubmult.mult.quantum_double.mult_poly_q_double`.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.from_expr"></a>

#### from\_expr

```python
def from_expr(expr)
```

Convert a polynomial to the quantum basis by repeatedly peeling off the leading monomial's
Schubert term; falls back to ``mul_expr`` on the identity if that fails.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.handle_sympoly"></a>

#### handle\_sympoly

```python
def handle_sympoly(other)
```

Evaluate symmetric-function coefficients to polynomials.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression: single ``x`` variables via ``mult_poly_q_double``, quantum factorial
elementary symmetric functions via ``elem_mul``, ``Add``/``Mul``/``Pow`` recursively, else as a coefficient.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QuantumDoubleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list or a polynomial expression.

<a id="schubmult.rings.schubert.quantum_double_schubert_ring.QDSx"></a>

#### QDSx

```python
def QDSx(x, genset=GeneratingSet("y"))
```

Construct a quantum double Schubert element in ``x`` with coefficient alphabet ``genset``
(a `GeneratingSet` or a label string); e.g. ``QDSx([3, 1, 2])`` is ``S^q_{312}(x; y)``.

<a id="schubmult.rings.schubert.quantum_schubert_ring"></a>

# schubmult.rings.schubert.quantum\_schubert\_ring

Quantum (single) Schubert polynomial ring: the ``QSx`` interface.

`QuantumSingleSchubertRing` is a `QuantumDoubleSchubertRing` with a zero
coefficient alphabet. This module also re-exports the quantum double and
parabolic quantum rings (``QDSx``, ``QPSx``, ``QPDSx``) for convenience.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing"></a>

## QuantumSingleSchubertRing Objects

```python
class QuantumSingleSchubertRing(QuantumDoubleSchubertRing)
```

The ring of quantum Schubert polynomials ``S^q_w(x)``; ``QSx`` is the standard instance.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.quantize"></a>

#### quantize

```python
def quantize(poly)
```

Quantize a polynomial: expand in classical Schubert polynomials, reinterpret each ``S_w`` as
the quantum ``S^q_w``, and expand back to a polynomial.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants: ``schubmult_q_fast`` when ``basis2`` is this ring, else ``schubmult_q_double_fast``.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product``.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply by an expression: single ``x`` variables via ``mult_poly_q``, ``Add``/``Mul``/``Pow``
recursively, anything else as a coefficient.

<a id="schubmult.rings.schubert.quantum_schubert_ring.QuantumSingleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list, a classical or parabolic element (converted
to the quantum basis), or a polynomial expression.

<a id="schubmult.rings.schubert.schubert_ring"></a>

# schubmult.rings.schubert.schubert\_ring

Ordinary (single) Schubert polynomial ring: the ``Sx`` interface.

`SingleSchubertRing` is a `DoubleSchubertRing` whose coefficient alphabet is
identically zero, so ``S_w(x; 0) = S_w(x)``. Products dispatch to the fast
integer kernel ``schubmult_py`` when both operands are single, and to
``schubmult_double`` when mixed with a genuinely double element.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing"></a>

## SingleSchubertRing Objects

```python
class SingleSchubertRing(DoubleSchubertRing)
```

The ring of ordinary Schubert polynomials ``S_w(x)``; ``Sx`` is the standard instance.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.cached_product"></a>

#### cached\_product

```python
@cache
def cached_product(u, v, basis2)
```

Structure constants of ``S_u * S_v``: integer ``schubmult_py`` when ``basis2`` is this ring,
else ``schubmult_double`` with ``y = 0``.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.cached_positive_product"></a>

#### cached\_positive\_product

```python
@cache
def cached_positive_product(u, v, basis2)
```

Same as ``cached_product`` (single coefficients are already nonnegative integers).

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.single_variable"></a>

#### single\_variable

```python
def single_variable(elem, varnum)
```

Multiply by ``x_varnum`` (non-equivariant Monk rule).

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.new"></a>

#### new

```python
def new(x)
```

Build an element from a permutation/Lehmer list or a polynomial expression.

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.elem_func"></a>

#### elem\_func

```python
@property
def elem_func()
```

`ElemSym` (non-factorial elementary symmetric function).

<a id="schubmult.rings.schubert.schubert_ring.SingleSchubertRing.divdiff"></a>

#### divdiff

```python
def divdiff(v, elem)
```

Apply the divided difference ``partial_v``: ``S_u -> S_{u v^{-1}}`` when length-additive, else 0.

<a id="schubmult.rings.schubert.separated_descents"></a>

# schubmult.rings.schubert.separated\_descents

Separated-descents ring: Schubert polynomials graded by an explicit number of variables.

`SeparatedDescentsRing` wraps a Schubert-family ring and indexes basis elements by
``(perm, num_vars)`` with ``num_vars >= max_descent(perm)``. The product of
``(u, p)`` and ``(v, q)`` places ``u`` in the first ``p`` variables and ``v`` in the
next ``q``, so descents of the two factors are separated. ``_sep_desc_mul`` computes
this (Samuel) as a twisted ordinary Schubert product by a dominant permutation. Also
carries experimental Pieri/coproduct routines (``pieri_formula``, ``coproduct_test``).

Not to be confused with `schubmult.mult.separated_descents`, which implements the
Fan-Guo-Xiong pipe-puzzle rule for double Grothendieck polynomials.

<a id="schubmult.rings.schubert.separated_descents.complete_sym_positional_perms_down"></a>

#### complete\_sym\_positional\_perms\_down

```python
def complete_sym_positional_perms_down(orig_perm, p, *k, hack_off=None)
```

Descent-side analogue of ``complete_sym_positional_perms``: all ``(perm, degree, sign)`` reachable
from ``orig_perm`` by up to ``p`` Bruhat *descents* swapping a fixed position in ``k`` (1-indexed)
with a still-untouched position. ``hack_off`` bounds the positions considered.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing"></a>

## SeparatedDescentsRing Objects

```python
class SeparatedDescentsRing(BaseSchubertRing)
```

Ring with basis ``(perm, num_vars)``; construct with ``SeparatedDescentsRing(Sx.ring)`` and call
as ``ring(perm, num_vars)``. See the module docstring.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.schub_ring"></a>

#### schub\_ring

```python
@property
def schub_ring()
```

The underlying Schubert-family ring used for products.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.pieri_formula"></a>

#### pieri\_formula

```python
def pieri_formula(p, elem)
```

Multiply ``elem`` by the single-row element ``(uncode([p]), 1)`` (adds one variable), via
``complete_sym_positional_perms_down``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(ring)
```

Wrap the Schubert-family ``ring`` (inherits its alphabets).

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.coproduct_test"></a>

#### coproduct\_test

```python
def coproduct_test(key)
```

Experimental coproduct of the basis element ``key = (perm, num_vars)``, computed by
triangular peeling of the leading code entry via ``pieri_formula``/``_single_coprod_test``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Ring product: scalars scale; otherwise each pair of basis elements multiplies via ``_sep_desc_mul``
and lands in degree ``deg1 + deg2`` (terms whose descents exceed that are dropped).

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

The ``SepDescSchubPoly`` display symbol for ``k = (perm, num_vars)``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRing.new"></a>

#### new

```python
def new(perm, deg=0)
```

Build ``(perm, deg)`` with ``deg`` raised to at least ``perm``'s last descent; a non-permutation
``perm`` is expanded in the underlying ring and each term given its minimal degree.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRingElement"></a>

## SeparatedDescentsRingElement Objects

```python
class SeparatedDescentsRingElement(BaseSchubertElement)
```

An element of a `SeparatedDescentsRing`: ``{(perm, num_vars): coeff}``.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRingElement.coproduct_test"></a>

#### coproduct\_test

```python
def coproduct_test()
```

Experimental coproduct; see `SeparatedDescentsRing.coproduct_test`.

<a id="schubmult.rings.schubert.separated_descents.SeparatedDescentsRingElement.as_ordered_terms"></a>

#### as\_ordered\_terms

```python
def as_ordered_terms(*_, **__)
```

Terms sorted by permutation length, permutation, then ``num_vars`` (sympy printing hook).

<a id="schubmult.rings.tensor_ring"></a>

# schubmult.rings.tensor\_ring

`TensorRing`: tensor products ``R_1 (x) ... (x) R_n`` of `BaseRing` instances.

Built with the ``@`` operator on rings (``Sx @ Sx``) or ``TensorRing(R1, R2, ...)``; nested
tensor rings are flattened. Keys are tuples ``(k_1, ..., k_n)`` of factor keys, multiplication
is factorwise, and the coproduct of a ring lands in ``R @ R``. Elements print as
``a # b``.

<a id="schubmult.rings.tensor_ring.TensorRing"></a>

## TensorRing Objects

```python
class TensorRing(BaseRing)
```

Tensor product of rings; keys are tuples of factor keys. See the module docstring.

<a id="schubmult.rings.tensor_ring.TensorRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Coproduct of the basis element ``k = (k_1, ..., k_n)`` into ``self @ self``.

Takes each factor's coproduct and interlaces them into flat keys
``(k_1^L, ..., k_n^L, k_1^R, ..., k_n^R)``.

<a id="schubmult.rings.tensor_ring.TensorRing.from_rc_graph_tensor"></a>

#### from\_rc\_graph\_tensor

```python
def from_rc_graph_tensor(rc_graph_tensor)
```

Pure tensor of the two factor rings' ``from_rc_graph`` images of a pair of RC graphs.

<a id="schubmult.rings.tensor_ring.TensorRing.__init__"></a>

#### \_\_init\_\_

```python
def __init__(*rings)
```

Tensor the given rings, flattening any that are themselves tensor rings.

<a id="schubmult.rings.tensor_ring.TensorRing.rings"></a>

#### rings

```python
@property
def rings()
```

The (flattened) tuple of tensor factors.

<a id="schubmult.rings.tensor_ring.TensorRing.rmul"></a>

#### rmul

```python
def rmul(elem1, elem2)
```

Scale every coefficient of ``elem1`` by the scalar ``elem2``.

<a id="schubmult.rings.tensor_ring.TensorRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Factorwise product: ``(a_1 (x) ... (x) a_n) * (b_1 (x) ... (x) b_n) = (a_1 b_1) (x) ... (x) (a_n b_n)``,
expanding each factor product in its own ring.

<a id="schubmult.rings.tensor_ring.TensorRing.cached_schubpoly"></a>

#### cached\_schubpoly

```python
@cache
def cached_schubpoly(k)
```

Product of the factor rings' polynomials for the key tuple ``k``.

<a id="schubmult.rings.tensor_ring.TensorRing.from_comp_ring"></a>

#### from\_comp\_ring

```python
def from_comp_ring(t)
```

Embed an element of one factor (or of a sub-tensor of factors) into this ring, filling the
other positions with their ``zero_monom`` (the identity).

<a id="schubmult.rings.tensor_ring.TensorRing.ext_multiply"></a>

#### ext\_multiply

```python
def ext_multiply(elem1, elem2)
```

External (tensor) product ``elem1 (x) elem2``: concatenates keys, flattening tensor-ring inputs.

<a id="schubmult.rings.tensor_ring.TensorRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(x)
```

A key tuple gives the corresponding basis element; anything else is parsed via ``from_expr``.

<a id="schubmult.rings.tensor_ring.TensorBasisElement"></a>

## TensorBasisElement Objects

```python
class TensorBasisElement(PrintingTerm)
```

Printing term for a tensor key; renders as ``a # b`` (str) or a tensor product (pretty/LaTeX).

<a id="schubmult.rings.tensor_ring.TensorRingElement"></a>

## TensorRingElement Objects

```python
class TensorRingElement(BaseRingElement)
```

Element of a `TensorRing`: a dict from key tuples to coefficients.

<a id="schubmult.rings.tensor_ring.TensorRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Coproduct into ``ring @ ring`` via `TensorRing.coproduct_on_basis`.

<a id="schubmult.rings.tensor_ring.TensorRingElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand to a commutative polynomial by multiplying out the factors' expansions (all factors
are assumed to live in disjoint or commuting variable sets).

<a id="schubmult.rings.thompson_algebra"></a>

# schubmult.rings.thompson\_algebra

`ThompsonAlgebra`: a noncommutative algebra on words in generators ``T_i`` and ``R_i``.

A monomial is a tuple of nonzero integers: ``i > 0`` stands for ``T_i`` and ``i < 0`` for
``R_{-i}``. Products are rewritten to a normal form by the rules in `ThompsonAlgebra._commute_pair`:
``T_i T_j = T_j T_{i+1}`` for ``i > j`` (the Thompson monoid relation), and ``T_i R_j`` moves
``R`` to the left with the index shifts (and one two-term case) given there.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebraElement"></a>

## ThompsonAlgebraElement Objects

```python
class ThompsonAlgebraElement(BaseRingElement)
```

Element of `ThompsonAlgebra`: a dict from normal-form words to coefficients.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra"></a>

## ThompsonAlgebra Objects

```python
class ThompsonAlgebra(BaseRing)
```

Algebra on words in ``T_i`` (positive index) and ``R_i`` (negative index). See the module docstring.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra.printing_term"></a>

#### printing\_term

```python
@cache
def printing_term(monomial)
```

Noncommutative product of the ``T_i``/``R_i`` symbols for the word.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra.new"></a>

#### new

```python
def new(x)
```

Build an element from a word (normalized via `_mul_monomials`), a number, or an existing element.

<a id="schubmult.rings.thompson_algebra.ThompsonAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Scalar multiplication, or the bilinear extension of `_mul_monomials`.

<a id="schubmult.symbolic"></a>

# schubmult.symbolic

Symbolic computation facade: fast SymEngine arithmetic with SymPy printing and polynomial domains.

Import ``Add``, ``Mul``, ``Pow``, ``S``, ``Symbol``, ``sympify``, ``expand`` from here rather than
from ``symengine``/``sympy`` directly; the SymPy versions are available under ``sympy_``-prefixed
names (``sympy_Add``, ``sympy_Mul``, ``sympify_sympy``, ``sympy_poly``). Also re-exports the
Schubert-polynomial helpers of `schubmult.symbolic.poly.schub_poly` and SymPy's
``EXRAW``/``CoercionFailed`` used by the ring domains. Generating sets live in
`schubmult.symbolic.poly.variables` (re-exported from `schubmult.symbolic.poly`).

<a id="schubmult.symbolic.common_polys"></a>

# schubmult.symbolic.common\_polys

Re-exports `schubmult.symbolic.poly.schub_poly` and the ``_vars``/``call_zvars``/``q_vector`` helpers.

<a id="schubmult.symbolic.functions"></a>

# schubmult.symbolic.functions

SymEngine-first wrappers (`expand`, `symbols`, `sympify`) that fall back to SymPy, plus small helpers.

<a id="schubmult.symbolic.functions.expand"></a>

#### expand

```python
def expand(obj, **kwargs)
```

Expand with SymEngine; use SymPy if keyword options are given or SymEngine fails.

<a id="schubmult.symbolic.functions.symbols"></a>

#### symbols

```python
def symbols(*args, **kwargs)
```

SymEngine ``symbols``.

<a id="schubmult.symbolic.functions.sympify"></a>

#### sympify

```python
def sympify(val)
```

SymEngine ``sympify``, falling back to SymPy for objects SymEngine cannot convert.

<a id="schubmult.symbolic.functions.is_of_func_type"></a>

#### is\_of\_func\_type

```python
def is_of_func_type(elem, typ)
```

``isinstance`` that also sees through SymEngine ``PyFunction`` wrappers around SymPy functions.

<a id="schubmult.symbolic.functions.expand_seq"></a>

#### expand\_seq

```python
def expand_seq(seq, genset)
```

The monomial ``genset[1]**seq[0] * genset[2]**seq[1] * ...`` (1-indexed generators).

<a id="schubmult.symbolic.functions.efficient_subs"></a>

#### efficient\_subs

```python
def efficient_subs(expr, subs_dict)
```

``expr.subs`` restricted to the entries of ``subs_dict`` that actually occur in ``expr``.

<a id="schubmult.symbolic.poly"></a>

# schubmult.symbolic.poly

Explicit polynomial machinery: generating sets (`variables`) and Schubert polynomial formulas (`schub_poly`).

<a id="schubmult.symbolic.poly.schub_poly"></a>

# schubmult.symbolic.poly.schub\_poly

Explicit symbolic formulas for (double) Schubert and Grothendieck polynomials.

The main entry points are `schubpoly` (double Schubert polynomial by the ``pull_out_var``
recursion), `schubpoly_from_elems` (Schubert polynomial as a sum of products of elementary
symmetric polynomials along theta-code v-paths, with a pluggable ``elem_func``),
`grothendieck_poly` (via isobaric divided differences), and the divided-difference operators
`div_diff`/`divide_out_diff`. ``_vars`` holds the default generating sets ``x``, ``y``, ``z``,
``q``. Everything here works on raw SymEngine/SymPy expressions; the ring classes in
`schubmult.rings` call these to expand basis elements.

<a id="schubmult.symbolic.poly.schub_poly.sv_posify"></a>

#### sv\_posify

```python
def sv_posify(val, var2)
```

Rewrite ``val`` in the differences ``var2[i+1] - var2[i]`` of consecutive variables (a
positivity-revealing form), by substituting ``var2[i] = var2[1] + r_1 + ... + r_{i-1}``,
simplifying, and mapping the ``r`` variables back.

<a id="schubmult.symbolic.poly.schub_poly.act"></a>

#### act

```python
def act(w, poly, genset)
```

Permute the variables of ``poly``: ``genset[i] -> genset[w(i)]``.

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_func"></a>

#### elem\_sym\_func

```python
def elem_sym_func(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2)
```

The double elementary symmetric factor attached to one step of the ``schubmult_double`` v-path
recursion: ``e_{k - udiff - vdiff}`` in the ``y`` variables fixed by ``u1 -> u2`` and the ``z``
variables selected by `call_zvars` for ``v1 -> v2``.

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_func_q"></a>

#### elem\_sym\_func\_q

```python
def elem_sym_func_q(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2)
```

Quantum-double variant of `elem_sym_func` (all ``k`` positions of ``u1``/``u2`` are compared).

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_poly_q"></a>

#### elem\_sym\_poly\_q

```python
def elem_sym_poly_q(p, k, varl1, varl2, q_var=_vars.q_var)
```

Quantum double elementary symmetric polynomial ``E_p^q(x_1..x_k; y)``: the usual recursion
plus the term ``q_{k-1} E_{p-2}(x_1..x_{k-2})``.

<a id="schubmult.symbolic.poly.schub_poly.complete_sym_poly"></a>

#### complete\_sym\_poly

```python
def complete_sym_poly(p, k, vrs, vrs2)
```

Factorial complete homogeneous symmetric polynomial ``h_p(vrs[0..k-1] | vrs2)``, computed by
splitting the variable set in half.

<a id="schubmult.symbolic.poly.schub_poly.elem_sym_poly"></a>

#### elem\_sym\_poly

```python
def elem_sym_poly(p, k, varl1, varl2, xstart=0, ystart=0)
```

Factorial elementary symmetric polynomial ``e_p(x_1 - y_1, ..., x_k - y_k)`` style sum over
``varl1[xstart:xstart+k]`` and ``varl2[ystart:]``, computed by a divide-and-conquer split of the
variables (the ``y`` offset shifts by the degree taken from the first half).

<a id="schubmult.symbolic.poly.schub_poly.call_zvars"></a>

#### call\_zvars

```python
@cache
def call_zvars(v1, v2, k, i, min_size=10)
```

Indices of the ``z`` variables entering the elementary symmetric factor for the v-path step
``v1 -> v2`` at position ``i`` with ``k`` variables (cached).

<a id="schubmult.symbolic.poly.schub_poly.q_vector"></a>

#### q\_vector

```python
def q_vector(q_exp, q_var=_vars.q_var)
```

Exponent vector of a monomial in the ``q`` variables (``q_1^a q_2^b -> [a, b]``); ``[]`` for 1,
``None`` if ``q_exp`` is not a ``q`` monomial.

<a id="schubmult.symbolic.poly.schub_poly.monom_sym"></a>

#### monom\_sym

```python
def monom_sym(partition, numvars, genset)
```

Monomial symmetric polynomial ``m_partition(genset[1..numvars])``.

<a id="schubmult.symbolic.poly.schub_poly.xreplace_genvars"></a>

#### xreplace\_genvars

```python
def xreplace_genvars(poly, vars1, vars2)
```

Replace the internal placeholder generating sets ``_vars.var_g1``/``var_g2`` with ``vars1``/``vars2``.

<a id="schubmult.symbolic.poly.schub_poly.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(poly, v1, v2)
```

The quotient ``(poly - poly|_{v1 -> v2}) / (v1 - v2)``, computed structurally on the expression
tree (so it is exact and needs no polynomial division). Objects may override via
``_eval_divide_out_diff``.

<a id="schubmult.symbolic.poly.schub_poly.split_up"></a>

#### split\_up

```python
def split_up(poly, v1, v2)
```

Write ``poly = a + (v1 - v2) * b`` with ``a = poly|_{v1 -> v2}``; returns ``(a, (v1 - v2, b))``.

<a id="schubmult.symbolic.poly.schub_poly.perm_act"></a>

#### perm\_act

```python
def perm_act(val, i, var2=None)
```

Swap ``var2[i]`` and ``var2[i+1]`` in ``val`` (the simple transposition ``s_i`` acting on variables).

<a id="schubmult.symbolic.poly.schub_poly.elem_func_func"></a>

#### elem\_func\_func

```python
def elem_func_func(k, i, v1, v2, vdiff, varl1, varl2, elem_func)
```

Single-sided version of `elem_sym_func` with a pluggable ``elem_func(p, k, xvars, zvars)``,
used by `schubpoly_from_elems`.

<a id="schubmult.symbolic.poly.schub_poly.elem_func_func_mul"></a>

#### elem\_func\_func\_mul

```python
def elem_func_func_mul(k, i, u1, u2, v1, v2, udiff, vdiff, varl1, varl2,
                       elem_func)
```

`elem_sym_func` with a pluggable ``elem_func`` in place of `elem_sym_poly`.

<a id="schubmult.symbolic.poly.schub_poly.schubpoly_from_elems"></a>

#### schubpoly\_from\_elems

```python
def schubpoly_from_elems(v, var_x=None, var_y=None, elem_func=None, mumu=None)
```

Schubert polynomial of ``v`` as a sum over v-paths of products of ``elem_func`` factors.

Uses the strict theta code of ``v^{-1}`` (or the code of the dominant ``mumu`` if given) and
the v-path dictionaries of `schubmult.utils.schub_lib.compute_vpathdicts`; each step
contributes ``elem_func(p, k, xvars, zvars)``. With ``elem_func = elem_sym_poly`` this is the
double Schubert polynomial; other choices give the SEM-basis expansion or, as in
`SchubertBasis.transition_word`, an encoding of the factors.

<a id="schubmult.symbolic.poly.schub_poly.schubpoly_classical_from_elems"></a>

#### schubpoly\_classical\_from\_elems

```python
def schubpoly_classical_from_elems(v, var_x=None, var_y=None, elem_func=None)
```

`schubpoly_from_elems` using the ordinary (non-strict) theta code of ``v^{-1}``.

<a id="schubmult.symbolic.poly.schub_poly.schubpoly"></a>

#### schubpoly

```python
def schubpoly(v, var2=None, var3=None, start_var=1)
```

Double Schubert polynomial ``S_v(var2; var3)`` by recursion on the last descent: pull out the
variable ``var2[n]`` (``n`` the last descent) via ``pull_out_var``, multiplying by factors
``(var2[n] - var3[p])``.

<a id="schubmult.symbolic.poly.schub_poly.div_diff"></a>

#### div\_diff

```python
def div_diff(poly, v1, v2)
```

Divided difference ``(poly - s(poly)) / (v1 - v2)`` where ``s`` swaps ``v1`` and ``v2``, computed
structurally on the expression tree. Objects may override via ``_eval_div_diff``.

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly_legacy"></a>

#### grothendieck\_poly\_legacy

```python
@cache
def grothendieck_poly_legacy(perm, x, y, beta, keep_as_schub=False)
```

Double Grothendieck polynomial by descending from ``w0`` (product of ``x (+) y`` factors) with
isobaric divided differences. Superseded by `grothendieck_poly`.

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly"></a>

#### grothendieck\_poly

```python
@cache
def grothendieck_poly(perm, x, y, beta, keep_as_schub=False)
```

Double Grothendieck polynomial ``G_perm(x; y)`` with parameter ``beta``, as an expression or
(``keep_as_schub``) as its double Schubert expansion. See `grothendieck_poly_with_ring`.

<a id="schubmult.symbolic.poly.schub_poly.dom_groth"></a>

#### dom\_groth

```python
@cache
def dom_groth(dom_perm, ring, beta)
```

Double Schubert expansion of the Grothendieck polynomial of a dominant permutation: builds
the product of factorial elementary symmetric factors row by row (from the code of
``dom_perm^{-1}``) with the ``1 + beta y`` twists.

<a id="schubmult.symbolic.poly.schub_poly.isobaric_strip_on_dschub_dict"></a>

#### isobaric\_strip\_on\_dschub\_dict

```python
def isobaric_strip_on_dschub_dict(start, length, perm_dict, coeff_genset,
                                  beta)
```

Apply one isobaric strip to a whole ``{perm: coeff}`` dict, folded.

Coefficients landing on the same permutation merge at every stage instead of
being carried per input basis element, mirroring ``compute_vpathdicts``.

<a id="schubmult.symbolic.poly.schub_poly.isobaric_strip_on_dschub"></a>

#### isobaric\_strip\_on\_dschub

```python
def isobaric_strip_on_dschub(start, length, schub_perm, ring, beta)
```

`isobaric_strip_on_dschub_dict` on a single basis element, returned as a ring element.

<a id="schubmult.symbolic.poly.schub_poly.apply_isobaric_to_schub_dict"></a>

#### apply\_isobaric\_to\_schub\_dict

```python
def apply_isobaric_to_schub_dict(diff_perm, perm_dict, coeff_genset, beta)
```

Fold every strip of ``diff_perm`` over the whole dict, merging between strips.

<a id="schubmult.symbolic.poly.schub_poly.apply_isobaric_to_schub"></a>

#### apply\_isobaric\_to\_schub

```python
@cache
def apply_isobaric_to_schub(diff_perm, schub_perm, ring, beta)
```

`apply_isobaric_to_schub_dict` on a single basis element, returned as a ring element.

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly_with_ring"></a>

#### grothendieck\_poly\_with\_ring

```python
@cache
def grothendieck_poly_with_ring(perm, ring, beta, keep_as_schub=False)
```

Double Grothendieck polynomial via the minimal dominant permutation above ``perm``: start
from `dom_groth` and apply the isobaric divided differences of ``perm^{-1} * dom_perm`` strip
by strip (`apply_isobaric_to_schub_dict`).

<a id="schubmult.symbolic.poly.schub_poly.grothendieck_poly2"></a>

#### grothendieck\_poly2

```python
@cache
def grothendieck_poly2(perm, x, y, beta, keep_as_schub=False)
```

Variant of `grothendieck_poly_legacy` with ``x - y - beta x y`` factors for ``w0``.

<a id="schubmult.symbolic.poly.schub_poly.to_groth"></a>

#### to\_groth

```python
def to_groth(val, x, y, beta)
```

Expand a polynomial in the double Grothendieck basis ``{perm: coeff}`` by triangular
elimination on monomials: peel off the lowest monomial ``x^c`` (lowest total degree, then lex),
subtract ``coeff * G_{uncode(c)}``, and recurse.

<a id="schubmult.symbolic.poly.schub_poly.to_groth_with_ring"></a>

#### to\_groth\_with\_ring

```python
def to_groth_with_ring(_val, ring, beta)
```

Expand a double Schubert ring element in the double Grothendieck basis.

Triangular elimination by length: for the smallest remaining permutation ``w``, apply the
isobaric divided differences of ``w`` and evaluate at ``x_i = -y_i / (1 + beta y_i)`` (the
point where all nontrivial Grothendieck polynomials vanish) to read off the coefficient of
``G_w``, then subtract ``coeff * G_w`` and repeat. Coefficients are simplified with SymPy.

<a id="schubmult.symbolic.poly.schub_poly.to_groth_with_ring_functional"></a>

#### to\_groth\_with\_ring\_functional

```python
def to_groth_with_ring_functional(_val, ring, beta)
```

`to_groth_with_ring` using the ring element's own ``isobaric_perm`` method.

<a id="schubmult.symbolic.poly.schub_poly.groth_dict_to_poly"></a>

#### groth\_dict\_to\_poly

```python
def groth_dict_to_poly(groth_dict, x, zz, beta)
```

Sum ``coeff * G_perm(x; zz)`` over a ``{perm: coeff}`` dict.

<a id="schubmult.symbolic.poly.schub_poly.schub_elem_to_groth_elem_dict"></a>

#### schub\_elem\_to\_groth\_elem\_dict

```python
@cache
def schub_elem_to_groth_elem_dict(the_perm, beta)
```

Signed count, by ``(inv, max_descent)``, of the permutations ``co_pipe_dream(rc).perm * w0`` over
RC graphs of ``the_perm``, weighted ``(-beta)^(inv difference)``: the Grothendieck-side image of a
Schubert basis element.

<a id="schubmult.symbolic.poly.schub_poly.schub_elem_sym_to_groth_elem_sym_dict"></a>

#### schub\_elem\_sym\_to\_groth\_elem\_sym\_dict

```python
@cache
def schub_elem_sym_to_groth_elem_sym_dict(p, k, beta)
```

`schub_elem_to_groth_elem_dict` for the Grassmannian permutation of ``e_p(x_1..x_k)``, i.e. the
expansion of the elementary symmetric polynomial into Grothendieck-Pieri pieces ``(inv, numvars)``.

<a id="schubmult.symbolic.poly.schub_poly.isobar_it"></a>

#### isobar\_it

```python
def isobar_it(i, genset, elem)
```

K-theoretic isobaric operator ``pi_i`` on a Schubert element: ``partial_i((1 + x_{i+1}) x_i * elem)``
via the nil-Hecke ring (``beta = 1``).

<a id="schubmult.symbolic.poly.schub_poly.lascoux_poly"></a>

#### lascoux\_poly

```python
def lascoux_poly(composition, genset)
```

Lascoux polynomial of a weak composition (``beta = 1``), expanded.

<a id="schubmult.symbolic.poly.schub_poly.groth_elem_as_schub_dict"></a>

#### groth\_elem\_as\_schub\_dict

```python
@cache
def groth_elem_as_schub_dict(perm, beta)
```

Schubert expansion ``{perm': coeff}`` of the Grothendieck polynomial ``G_perm`` (via
``WCGraph.groth_to_schub``).

<a id="schubmult.symbolic.poly.schub_poly.groth_mul_full"></a>

#### groth\_mul\_full

```python
def groth_mul_full(perm_dict, p2, _x, _zz, beta)
```

Multiply a Grothendieck expansion ``perm_dict`` by ``G_p2``: expand ``G_p2`` in Schubert
polynomials and push each through `schub_dict_to_groth_dict`.

<a id="schubmult.symbolic.poly.schub_poly.groth_mul_full_with_ring"></a>

#### groth\_mul\_full\_with\_ring

```python
def groth_mul_full_with_ring(perm_dict, p2, ring, beta)
```

`groth_mul_full` using the ring-aware `schub_dict_to_groth_dict_with_ring`.

<a id="schubmult.symbolic.poly.schub_poly.schub_dict_to_groth_dict"></a>

#### schub\_dict\_to\_groth\_dict

```python
def schub_dict_to_groth_dict(base_groth, schub_dict, beta)
```

Multiply the Grothendieck expansion ``base_groth`` by the Schubert polynomial with expansion
``schub_dict``, returning a Grothendieck expansion.

Writes the Schubert polynomial in the CEM (elementary symmetric) basis, converts each
``e_p(x_1..x_k)`` factor to Grothendieck-Pieri pieces with
`schub_elem_sym_to_groth_elem_sym_dict`, and applies ``groth_pieri_mul`` factor by factor.

<a id="schubmult.symbolic.poly.schub_poly.schub_dict_to_groth_dict_with_ring"></a>

#### schub\_dict\_to\_groth\_dict\_with\_ring

```python
def schub_dict_to_groth_dict_with_ring(base_groth, schub_dict, ring, beta)
```

`schub_dict_to_groth_dict` for a specific ``ring`` (uses ``ring.in_CEM_basis`` and
``ring.is_elem_mul_type`` to recognize elementary symmetric factors).

<a id="schubmult.symbolic.poly.variables"></a>

# schubmult.symbolic.poly.variables

Generating sets: indexed families of variables ``x_0, x_1, x_2, ...`` used as ring generators.

`GeneratingSet("x")` interns ``DEF_GENSET_SIZE`` symbols ``x_0..x_99``; ``gs[i]`` is the symbol
``x_i``, and the polynomial variables are ``x_1, x_2, ...`` (``x_0`` is unused), so an exponent
tuple ``(a_1, ..., a_n)`` means ``x_1^{a_1} ... x_n^{a_n}``. `MaskedGeneratingSet` hides a set of indices of a base set,
`CustomGeneratingSet` wraps an arbitrary sequence of expressions, and `ZeroGeneratingSet`
returns 0 for every index (used for single Schubert polynomials as a degenerate coefficient
set). `genset_dict_from_expr` converts a polynomial expression into ``{exponent_tuple: coeff}``.

<a id="schubmult.symbolic.poly.variables.GeneratingSet_base"></a>

## GeneratingSet\_base Objects

```python
class GeneratingSet_base()
```

Interface for generating sets: indexing, length, ``index(symbol)`` (``-1`` if absent), and ``label``.

<a id="schubmult.symbolic.poly.variables.ZeroGeneratingSet"></a>

## ZeroGeneratingSet Objects

```python
class ZeroGeneratingSet(GeneratingSet_base)
```

A generating set every entry of which is ``0``; contains no symbols.

<a id="schubmult.symbolic.poly.variables.GeneratingSet"></a>

## GeneratingSet Objects

```python
class GeneratingSet(GeneratingSet_base)
```

The interned family ``name_0, name_1, ...``; ``gs[i]`` is the symbol ``name_i`` and ``gs(i)`` is ``gs[i - 1]``.

<a id="schubmult.symbolic.poly.variables.GeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.GeneratingSet.label"></a>

#### label

```python
@property
def label()
```

The variable name, e.g. ``"x"``.

<a id="schubmult.symbolic.poly.variables.GeneratingSet.index"></a>

#### index

```python
def index(v)
```

Position of the symbol ``v`` in this set, or ``-1``.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet"></a>

## MaskedGeneratingSet Objects

```python
class MaskedGeneratingSet(GeneratingSet_base)
```

A base generating set with the (1-indexed) positions in ``index_mask`` removed and the rest
renumbered consecutively; ``complement()`` gives the set of the masked variables instead.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.base_genset"></a>

#### base\_genset

```python
@property
def base_genset()
```

The underlying unmasked generating set.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.index_mask"></a>

#### index\_mask

```python
@property
def index_mask()
```

Sorted tuple of the hidden 1-indexed positions.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.complement"></a>

#### complement

```python
def complement()
```

The masked set on the complementary positions.

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.CustomGeneratingSet"></a>

## CustomGeneratingSet Objects

```python
class CustomGeneratingSet(GeneratingSet_base)
```

A generating set over an explicit sequence of expressions (sympified on construction).

<a id="schubmult.symbolic.poly.variables.CustomGeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.NotEnoughGeneratorsError"></a>

## NotEnoughGeneratorsError Objects

```python
class NotEnoughGeneratorsError(ValueError)
```

Raised when an operation needs more generators than a generating set provides.

<a id="schubmult.symbolic.poly.variables.poly_genset"></a>

#### poly\_genset

```python
@cache
def poly_genset(v: str)
```

``GeneratingSet(v)``, or a `ZeroGeneratingSet` for the sentinels ``ZeroVar``/``NoneVar``.

<a id="schubmult.symbolic.poly.variables.genset_dict_from_expr"></a>

#### genset\_dict\_from\_expr

```python
def genset_dict_from_expr(expr, genset, length=None)
```

Write a polynomial in the generators of ``genset`` as ``{exponent_tuple: coeff}``.

Exponent tuples are 0-indexed by generator position ``genset(i) -> tuple[i - 1]`` and have
length ``length`` (default: the largest generator index present). Factors free of the
generators go into the coefficient; a factor mixing generators with other symbols raises.

<a id="schubmult.symbolic.symmetric_polynomials"></a>

# schubmult.symbolic.symmetric\_polynomials

Symbolic (factorial) elementary and complete symmetric polynomials as SymPy functions.

``E(p, k, *vars)``/``e`` and ``H(p, k, *vars)``/``h`` are the elementary and complete symmetric
polynomials of degree ``p`` in the first ``k`` generators, kept unevaluated so Schubert
expansions can be written in the SEM basis; the ``Factorial*`` variants carry a second
variable set. `functions` holds canonicalization and variable-splitting utilities, and
`qelem_sym` the quantum elementary symmetric polynomials.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym"></a>

# schubmult.symbolic.symmetric\_polynomials.complete\_sym

Unevaluated (factorial) complete homogeneous symmetric polynomials as SymPy function atoms.

``H(p, k, xvars, yvars)`` (alias `FactorialCompleteSym`) is the factorial complete symmetric
polynomial of degree ``p`` in ``k`` generators with ``p + k - 1`` coefficient variables;
``h(p, k, xvars)`` (alias `CompleteSym`) is the ordinary ``h_p(x_1..x_k)``. The two families
are related by the duality ``H(p, k; x, y) = (-1)^p E(p, k + 1 - p; y, x)`` (`H.to_elem_sym`,
`H.from_elem_sym`), and divided differences are computed by passing through `E`.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.CompleteSym_base"></a>

## CompleteSym\_base Objects

```python
class CompleteSym_base(Function)
```

Common behavior for `H` and `h`; ``expand_func`` evaluates via `complete_sym_poly`.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H"></a>

## H Objects

```python
class H(CompleteSym_base)
```

Factorial complete symmetric polynomial ``H(p, k, xvars, yvars)``; see the module docstring.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.to_elem_sym"></a>

#### to\_elem\_sym

```python
def to_elem_sym()
```

``(-1)^p E(p, k + 1 - p; yvars, xvars)``: the same polynomial as a factorial elementary symmetric atom.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.from_elem_sym"></a>

#### from\_elem\_sym

```python
@classmethod
def from_elem_sym(cls, elem, sign=False)
```

Inverse of `to_elem_sym`: the ``H`` atom equal to the `E` atom ``elem`` (with the ``(-1)^p`` if ``sign``).

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(vars1, vars2=None)
```

``H_p(all) = sum_i H_i(vars1) H_{p-i}(rest)`` with the coefficient variables split accordingly.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(v1, v2)
```

`E.divide_out_diff` transported through `to_elem_sym`/`from_elem_sym`.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.from_expr_elem_sym"></a>

#### from\_expr\_elem\_sym

```python
@staticmethod
def from_expr_elem_sym(expr)
```

Replace every `E` atom in ``expr`` by the equal `H` atom.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.to_expr_elem_sym"></a>

#### to\_expr\_elem\_sym

```python
@staticmethod
def to_expr_elem_sym(expr)
```

Replace every `H` atom in ``expr`` by the equal `E` atom.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.H.div_diff"></a>

#### div\_diff

```python
def div_diff(v1, v2)
```

Divided difference, computed on the `E` side and converted back.

<a id="schubmult.symbolic.symmetric_polynomials.complete_sym.h"></a>

## h Objects

```python
class h(CompleteSym_base)
```

Ordinary complete homogeneous symmetric polynomial ``h(p, k, xvars)``.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym"></a>

# schubmult.symbolic.symmetric\_polynomials.elem\_sym

Unevaluated (factorial) elementary symmetric polynomials as SymPy function atoms.

``E(p, k, xvars, yvars)`` (alias `FactorialElemSym`) is the factorial elementary symmetric
polynomial of degree ``p`` in the ``k`` generators ``xvars`` with coefficient variables
``yvars`` (``k + 1 - p`` of them are used); ``e(p, k, xvars)`` (alias `ElemSym`) is the ordinary
``e_p(x_1..x_k)``. Both stay symbolic so Schubert polynomials can be manipulated in the SEM/CEM
bases; ``expand_func`` evaluates them via `schubmult.symbolic.poly.schub_poly.elem_sym_poly`.
They implement divided differences (`div_diff`, `divide_out_diff`), variable splitting
(`split_out_vars`), and canonicalize on construction (``E(p, k, ...) = 0`` if ``p > k``, ``1`` if
``p == 0``, and a shared variable between the two sets cancels).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base"></a>

## ElemSym\_base Objects

```python
class ElemSym_base(Function)
```

Common behavior for `E` and `e`: substitution acts only on the variable arguments, and the
``degree``/``numvars``/``genvars``/``coeffvars`` accessors expose the parameters.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.degree"></a>

#### degree

```python
@property
def degree()
```

``p``.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.numvars"></a>

#### numvars

```python
@property
def numvars()
```

``k``, the number of generators.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.genvars"></a>

#### genvars

```python
@property
def genvars()
```

The ``x`` variables (sorted tuple).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.ElemSym_base.coeffvars"></a>

#### coeffvars

```python
@property
def coeffvars()
```

The ``y`` (coefficient) variables.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E"></a>

## E Objects

```python
class E(ElemSym_base)
```

Factorial elementary symmetric polynomial ``E(p, k, xvars, yvars)``; see the module docstring.

Variables may be passed as two iterables or flattened (``k`` x's followed by ``k + 1 - p`` y's).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.cauchy"></a>

#### cauchy

```python
@staticmethod
def cauchy(fnc, genset)
```

Rewrite ``fnc`` so its coefficient variables are the initial segment of ``genset``, using
``E(p,k;..y..) = E(p,k;..y'..) + (y - y') E(p-1,k-1;..y'..)`` one variable at a time.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(vars1, vars2=None)
```

Split the generators into ``vars1`` and the rest: ``E(p, k) = sum_i E(i, k1; ..) E(p - i, k2; ..)``
with the coefficient variables distributed accordingly. With ``vars1=None`` the split is
made on the coefficient variables ``vars2`` instead.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(v1, v2)
```

``(self - self|_{v1 -> v2}) / (v1 - v2)`` in closed form: removing a generator lowers ``p`` and
``k`` by one; a coefficient variable gives the corresponding signed term.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.div_diff"></a>

#### div\_diff

```python
def div_diff(v1, v2)
```

Divided difference ``partial_{v1, v2}`` in closed form (antisymmetric in ``v1``, ``v2``).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.E.pull_out_vars"></a>

#### pull\_out\_vars

```python
def pull_out_vars(var1, var2, min_degree=1)
```

Write ``self = self|_{var1 -> var2} + (var1 - var2) * divide_out_diff(var1, var2)`` when ``var1``
is a generator and ``var2`` a coefficient variable (and ``p >= min_degree``).

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e"></a>

## e Objects

```python
class e(ElemSym_base)
```

Ordinary elementary symmetric polynomial ``e(p, k, xvars)``; see the module docstring.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(vars1, vars2=None)
```

``e_p(all) = sum_i e_i(vars1) e_{p-i}(rest)``.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.coeffvars"></a>

#### coeffvars

```python
@property
def coeffvars()
```

No coefficient variables: a `ZeroGeneratingSet`.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.divide_out_diff"></a>

#### divide\_out\_diff

```python
def divide_out_diff(v1, v2)
```

``e_{p-1}`` of the remaining generators if ``v1`` is a generator, else 0.

<a id="schubmult.symbolic.symmetric_polynomials.elem_sym.e.div_diff"></a>

#### div\_diff

```python
def div_diff(v1, v2)
```

Divided difference: ``+/- e_{p-1}`` of the remaining generators, or 0 if neither variable is a generator.

<a id="schubmult.symbolic.symmetric_polynomials.functions"></a>

# schubmult.symbolic.symmetric\_polynomials.functions

Expression-level utilities for symbolic elementary symmetric polynomials.

The accessors `genvars`/`coeffvars`/`degree`/`numvars` see through SymEngine ``PyFunction``
wrappers; `split_out_vars`/`pull_out_vars` map the corresponding `E` methods over a whole
expression tree; `canonicalize_elem_syms` rewrites products of ``E`` factors into a normal form
(each factor of full degree ``p == k``, grouped by first coefficient variable).

<a id="schubmult.symbolic.symmetric_polynomials.functions.genvars"></a>

#### genvars

```python
def genvars(obj)
```

``obj.genvars``, unwrapping a SymEngine ``PyFunction`` if needed.

<a id="schubmult.symbolic.symmetric_polynomials.functions.coeffvars"></a>

#### coeffvars

```python
def coeffvars(obj)
```

``obj.coeffvars``, unwrapping a SymEngine ``PyFunction`` if needed.

<a id="schubmult.symbolic.symmetric_polynomials.functions.degree"></a>

#### degree

```python
def degree(obj)
```

The degree ``p`` of an elementary symmetric atom (unwrapping if needed).

<a id="schubmult.symbolic.symmetric_polynomials.functions.numvars"></a>

#### numvars

```python
def numvars(obj)
```

The variable count ``k`` of an elementary symmetric atom (unwrapping if needed).

<a id="schubmult.symbolic.symmetric_polynomials.functions.canonicalize_elem_syms"></a>

#### canonicalize\_elem\_syms

```python
def canonicalize_elem_syms(expr, combine_equal=False)
```

Normal form for expressions in `FactorialElemSym`: split every factor with ``p < k`` in half
until all factors have ``p == k``, then within each product regroup factors sharing a first
coefficient variable (merging them into one factor if ``combine_equal``).

<a id="schubmult.symbolic.symmetric_polynomials.functions.canonicalize_elem_syms_coeff"></a>

#### canonicalize\_elem\_syms\_coeff

```python
def canonicalize_elem_syms_coeff(expr, combine_equal=False)
```

`canonicalize_elem_syms` splitting on coefficient variables instead of generators.

<a id="schubmult.symbolic.symmetric_polynomials.functions.split_out_vars"></a>

#### split\_out\_vars

```python
def split_out_vars(expr, vars1, vars2)
```

Apply ``split_out_vars(vars1, vars2)`` to every elementary symmetric atom in ``expr``.

<a id="schubmult.symbolic.symmetric_polynomials.functions.pull_out_vars"></a>

#### pull\_out\_vars

```python
def pull_out_vars(expr, var1, var2, min_degree=1)
```

Apply ``pull_out_vars(var1, var2, min_degree)`` to every elementary symmetric atom in ``expr``.

<a id="schubmult.symbolic.symmetric_polynomials.functions.elem_sym_unify"></a>

#### elem\_sym\_unify

```python
def elem_sym_unify(expr, arg=None)
```

Recursively walk ``expr`` (currently a structural no-op; the pattern-based unification is
commented out).

<a id="schubmult.symbolic.symmetric_polynomials.qelem_sym"></a>

# schubmult.symbolic.symmetric\_polynomials.qelem\_sym

Quantum factorial elementary symmetric polynomials ``E_q(p, k, xvars, yvars)`` as SymPy atoms.

The quantum deformation adds ``q_i`` terms for adjacent pairs ``x_i, x_{i+1}`` of the *positional*
generators (positions taken from ``x_var``), so the variable set need not be an initial segment.
``expand_func`` evaluates via `elem_sym_positional_poly_q`. Alias: `QFactorialElemSym`.

<a id="schubmult.symbolic.symmetric_polynomials.qelem_sym.E_q"></a>

## E\_q Objects

```python
class E_q(ElemSym_base)
```

Quantum factorial elementary symmetric atom; same canonicalization as `E`, plus ``x_var``/``q_var``
generating sets fixing the positions of the generators and the ``q`` parameters.

<a id="schubmult.symbolic.symmetric_polynomials.qelem_sym.elem_sym_positional_poly_q"></a>

#### elem\_sym\_positional\_poly\_q

```python
def elem_sym_positional_poly_q(p,
                               k,
                               varl1,
                               varl2,
                               x_var=GeneratingSet("x"),
                               q_var=GeneratingSet("q"))
```

Quantum factorial elementary symmetric polynomial of degree ``p`` in the generators ``varl1[:k]``.

Recursion on the last generator ``x_l``: the classical two terms plus, when ``x_{l-1}`` (resp.
``x_{l+1}``) is also among the generators, ``q_{l-1}`` (resp. ``q_l``) times the degree ``p - 2``
polynomial with that pair removed.

<a id="schubmult.utils"></a>

# schubmult.utils

Support utilities: permutation/list helpers (`perm_utils`), the multiplication kernels' combinatorial
core (`schub_lib`), CLI argument parsing (`argparse`), coefficient parsing, logging, and grid printing.

<a id="schubmult.utils._grid_print"></a>

# schubmult.utils.\_grid\_print

`GridPrint`: a SymPy ``Printable`` mixin that renders a 2-D grid (``rows``, ``cols``, ``self[i, j]``)
as an aligned table for str/pretty/LaTeX output; used by RC graphs, BPDs and similar diagrams.

<a id="schubmult.utils._grid_print.GridPrint"></a>

## GridPrint Objects

```python
class GridPrint(Printable)
```

Mixin: subclasses provide ``rows``, ``cols``, ``__getitem__((i, j))`` and ``_display_name``.

<a id="schubmult.utils._mul_utils"></a>

# schubmult.utils.\_mul\_utils

Dict-level helpers for ring multiplication and tensor products of ``{key: coeff}`` expansions.

<a id="schubmult.utils.argparse"></a>

# schubmult.utils.argparse

Shared command-line parsing for the ``schubmult_*`` scripts (`schub_argparse`).

<a id="schubmult.utils.argparse.schub_argparse"></a>

#### schub\_argparse

```python
def schub_argparse(prog_name,
                   description,
                   argv,
                   quantum=False,
                   yz=False,
                   coprod=True)
```

Parse the common CLI of the ``schubmult_*`` scripts and return ``(args, formatter)``.

Permutations are given as space-separated integers separated by ``-`` (or Lehmer codes with
``--code``); ``quantum`` adds the ``--parabolic`` options, ``yz`` the double-variable and
``--display-positive`` options, ``coprod`` the ``--coprod`` mode. ``--display-mode`` selects the
``formatter`` (a callable rendering expressions as LaTeX, pretty, basic, sympy, or ``None`` for
raw). Hidden ``-g`` dumps the parsed arguments to a JSON file for the script test suite and
exits. Also initializes SymPy printing and logging.

<a id="schubmult.utils.logging"></a>

# schubmult.utils.logging

Thin wrappers around the standard ``logging`` module.

<a id="schubmult.utils.logging.init_logging"></a>

#### init\_logging

```python
def init_logging(debug=False)
```

Configure root logging at DEBUG (if ``debug``) or ERROR with a timestamped ``file:line`` format.

<a id="schubmult.utils.logging.get_logger"></a>

#### get\_logger

```python
def get_logger(name)
```

``logging.getLogger(name)``.

<a id="schubmult.utils.parsing"></a>

# schubmult.utils.parsing

Parsing of coefficient expressions from the command line.

<a id="schubmult.utils.parsing.parse_coeff"></a>

#### parse\_coeff

```python
def parse_coeff(coeff_str, latex=False)
```

Sympify ``coeff_str`` and map any symbol ``name_i`` to ``GeneratingSet(name)[i]`` so the result
uses the package's interned variables. LaTeX input is not yet supported (returns ``None``).

<a id="schubmult.utils.perm_utils"></a>

# schubmult.utils.perm\_utils

Low-level helpers on permutations given as lists/tuples, plus composition and reduced-word utilities.

Most functions here are used by the multiplication kernels in `schubmult.mult` and by
`schubmult.combinatorics.permutation.Permutation`; `add_perm_dict` is the standard way to merge
``{key: coeff}`` expansions.

<a id="schubmult.utils.perm_utils.permtrim_list"></a>

#### permtrim\_list

```python
def permtrim_list(perm)
```

Strip trailing fixed points ``perm[L-1] == L`` from a list in place and return it.

<a id="schubmult.utils.perm_utils.has_bruhat_descent"></a>

#### has\_bruhat\_descent

```python
def has_bruhat_descent(perm, i, j)
```

Check if perm has a Bruhat descent from position i to j.

Optimized version assuming perm is a Permutation object with direct indexing.

<a id="schubmult.utils.perm_utils.count_bruhat"></a>

#### count\_bruhat

```python
def count_bruhat(perm, i, j)
```

Signed length change ``inv(perm * t_{ij}) - inv(perm)`` for the transposition of positions ``i < j``.

<a id="schubmult.utils.perm_utils.has_bruhat_ascent"></a>

#### has\_bruhat\_ascent

```python
def has_bruhat_ascent(perm, i, j)
```

Check if perm has a Bruhat ascent from position i to j.

Optimized version assuming perm is a Permutation object with direct indexing.

<a id="schubmult.utils.perm_utils.omega"></a>

#### omega

```python
def omega(i, qv)
```

``i``-th entry (1-indexed) of the Cartan-matrix image of the q-exponent vector ``qv``:
``2 qv[i] - qv[i-1] - qv[i+1]`` with boundary conventions. Used to convert quantum ``q``
monomials to weights in the parabolic quantum product.

<a id="schubmult.utils.perm_utils.sg"></a>

#### sg

```python
def sg(i, w)
```

1 if ``w`` has a descent at 0-indexed position ``i``, else 0.

<a id="schubmult.utils.perm_utils.count_less_than"></a>

#### count\_less\_than

```python
def count_less_than(arr, val)
```

Number of leading entries of the sorted list ``arr`` that are ``< val``.

<a id="schubmult.utils.perm_utils.artin_sequences"></a>

#### artin\_sequences

```python
def artin_sequences(n)
```

All tuples ``(a_1, ..., a_n)`` with ``0 <= a_i <= n + 1 - i`` (Lehmer codes of ``S_{n+1}``).

<a id="schubmult.utils.perm_utils.weak_compositions"></a>

#### weak\_compositions

```python
def weak_compositions(length, max_degree)
```

All tuples of the given ``length`` with entries in ``0..max_degree``.

<a id="schubmult.utils.perm_utils.is_parabolic"></a>

#### is\_parabolic

```python
def is_parabolic(w, parabolic_index)
```

Whether ``w`` has no descent at any of the (1-indexed) positions in ``parabolic_index``.

<a id="schubmult.utils.perm_utils.add_perm_dict"></a>

#### add\_perm\_dict

```python
def add_perm_dict(d1, d2)
```

Return ``d1 + d2`` as coefficient dicts (keys merged, values added).

<a id="schubmult.utils.perm_utils.add_perm_dict_with_coeff"></a>

#### add\_perm\_dict\_with\_coeff

```python
def add_perm_dict_with_coeff(d1, d2, coeff)
```

Return ``d1 + coeff * d2`` as coefficient dicts.

<a id="schubmult.utils.perm_utils.p_trans"></a>

#### p\_trans

```python
def p_trans(part)
```

Conjugate (transpose) of a partition given as a weakly decreasing list; ``[0]`` for the empty partition.

<a id="schubmult.utils.perm_utils.mu_A"></a>

#### mu\_A

```python
def mu_A(mu, A)
```

The partition whose conjugate consists of the columns of ``mu`` indexed by ``A`` (0-indexed).

<a id="schubmult.utils.perm_utils.get_cycles"></a>

#### get\_cycles

```python
def get_cycles(perm)
```

``perm.get_cycles()``.

<a id="schubmult.utils.perm_utils.old_code"></a>

#### old\_code

```python
def old_code(perm)
```

Lehmer code of a permutation list computed by successive deletion from ``[1..L]``.

<a id="schubmult.utils.perm_utils.cyclic_sort"></a>

#### cyclic\_sort

```python
def cyclic_sort(L)
```

Rotate the list so its maximum is last.

<a id="schubmult.utils.perm_utils.cyclic_sort_min"></a>

#### cyclic\_sort\_min

```python
def cyclic_sort_min(L)
```

Rotate the list so its minimum is first.

<a id="schubmult.utils.perm_utils.h_vector"></a>

#### h\_vector

```python
def h_vector(q_vector)
```

Positions (1-indexed) where the vector strictly increases, up to its first decrease.

<a id="schubmult.utils.perm_utils.l_vector"></a>

#### l\_vector

```python
def l_vector(q_vector)
```

Find l_j = last position where d equals j (where d decreases from j to j-1).

<a id="schubmult.utils.perm_utils.tau_d"></a>

#### tau\_d

```python
def tau_d(d)
```

Partial permutation built from `h_vector`/`l_vector` of ``d`` (``tau[l_i - i] = h_i``), completed by
``Permutation.from_partial``.

<a id="schubmult.utils.perm_utils.phi_d"></a>

#### phi\_d

```python
def phi_d(d)
```

Companion of `tau_d` with the shifted placement ``phi[l_i - 1 - i] = h_i``.

<a id="schubmult.utils.perm_utils.conjugate_weak_composition"></a>

#### conjugate\_weak\_composition

```python
def conjugate_weak_composition(comp)
```

Compute the conjugate of a weak composition.

The conjugate of a weak composition α = (α₁, α₂, ..., αₙ) is the weak composition
β where βⱼ = |{i : αᵢ ≥ j}|, i.e., βⱼ counts how many parts of α are at least j.

This is equivalent to transposing the Ferrers diagram of the composition.

**Arguments**:

- `comp` - A sequence (list, tuple) of non-negative integers representing a weak composition.
  

**Returns**:

  A tuple representing the conjugate weak composition.
  

**Examples**:

  >>> conjugate_weak_composition([3, 1, 0, 2])
  (3, 2, 1)
  >>> conjugate_weak_composition([4, 2, 1])
  (3, 2, 1, 1)
  >>> conjugate_weak_composition([])
  ()
  >>> conjugate_weak_composition([0, 0, 0])
  ()

<a id="schubmult.utils.perm_utils.find_reduced_fail"></a>

#### find\_reduced\_fail

```python
def find_reduced_fail(word, inserted)
```

After changing letter ``inserted`` of a word, find the other position carrying the same root
(the letter whose deletion would make the word reduced again), or ``None``.

<a id="schubmult.utils.perm_utils.is_reduced"></a>

#### is\_reduced

```python
def is_reduced(word)
```

Whether the word of simple reflections is reduced (``inv`` of its product equals its length).

<a id="schubmult.utils.perm_utils.little_bump_pos"></a>

#### little\_bump\_pos

```python
def little_bump_pos(word, index)
```

Little bump at position ``index``: decrement that letter (increment if it is 1), and while the
word is not reduced, repeat at the letter found by `find_reduced_fail`.

<a id="schubmult.utils.perm_utils.little_bump"></a>

#### little\_bump

```python
def little_bump(word, i, j)
```

Little bump of a reduced word at the letter whose right root is the inversion ``(i, j)``.

<a id="schubmult.utils.perm_utils.little_zero"></a>

#### little\_zero

```python
def little_zero(word, length)
```

Repeatedly Little-bump at the last descent until the product's code has fewer than ``length``
entries (Little's map toward a smaller permutation).

<a id="schubmult.utils.schub_lib"></a>

# schubmult.utils.schub\_lib

Combinatorial kernels behind the Schubert multiplication algorithms.

The central objects are the *v-path dictionaries* (`compute_vpathdicts`): for a theta code
``th`` and target permutation ``vmu`` they record, level by level, the Bruhat-descent moves
(`kdown_perms`) along which the Schubert product is accumulated in `schubmult.mult.single` and
its double/quantum relatives. Also here: the ``elem_sym_perms*``/``complete_sym_perms*`` families
(Pieri-type moves for multiplying by elementary/complete symmetric polynomials, including quantum
and Grothendieck variants), `pull_out_var` (the recursion behind `schubpoly`), coefficient
reductions (`reduce_coeff`, `reduce_descents`, ``try_reduce_*``) that shrink a triple ``(u, v, w)``
before computing an LR coefficient, and the parabolic ``q``-vector checks.

<a id="schubmult.utils.schub_lib.double_elem_sym_q"></a>

#### double\_elem\_sym\_q

```python
def double_elem_sym_q(u, p1, p2, k, q_var=q_var)
```

Pairs of consecutive quantum Pieri moves ``u -> perm1 -> perm2`` (degrees ``p1`` then ``p2`` in ``k``
variables) whose cycles are compatible; returned as ``{(perm1, udiff1, q1): [(perm2, udiff2, q2), ...]}``.

<a id="schubmult.utils.schub_lib.will_formula_work"></a>

#### will\_formula\_work

```python
def will_formula_work(u, v)
```

Whether the fast dominant-descent formula applies to the pair ``(u, v)``: ``v^{-1} * mu_v`` must
already be the identity (no descents to reduce).

<a id="schubmult.utils.schub_lib.try_reduce_u"></a>

#### try\_reduce\_u

```python
def try_reduce_u(u, v, w)
```

Move a zero in the code of ``u`` past a nonzero entry by swapping adjacent positions in ``(u, w)``
or ``(u, v)`` where the LR coefficient is unchanged, aiming for ``u`` to one-dominate ``w``.

<a id="schubmult.utils.schub_lib.reduce_descents"></a>

#### reduce\_descents

```python
def reduce_descents(u, v, w)
```

Cancel common descents of ``w`` with ``v`` (or ``u``) by simultaneous adjacent swaps, stopping once
one of the fast-path conditions holds.

<a id="schubmult.utils.schub_lib.is_reducible"></a>

#### is\_reducible

```python
def is_reducible(v)
```

Whether the code of ``v`` has no nonzero entry after a zero (so ``v`` is a shifted dominant-like shape).

<a id="schubmult.utils.schub_lib.try_reduce_v"></a>

#### try\_reduce\_v

```python
def try_reduce_v(u, v, w)
```

`try_reduce_u` with the roles of ``u`` and ``v`` exchanged, aiming for `is_reducible`.

<a id="schubmult.utils.schub_lib.reduce_coeff"></a>

#### reduce\_coeff

```python
def reduce_coeff(u, v, w)
```

Reduce the LR triple ``(u, v, w)`` by dividing out the dominant permutations of the theta codes of
``u^{-1}`` and ``v^{-1}``: if ``w`` decomposes compatibly, return the smaller triple
``(u mu_A, v mu_B, w')`` with the same coefficient; otherwise return the input.

<a id="schubmult.utils.schub_lib.pull_out_var"></a>

#### pull\_out\_var

```python
def pull_out_var(vnum, v)
```

All ways to factor the variable ``x_vnum`` out of ``S_v``: returns pairs ``(indices, v')`` such that
``S_v = sum prod_{p in indices} (x_vnum - y_p) * S_{v'}`` with ``x_vnum`` absent from ``S_{v'}``.

<a id="schubmult.utils.schub_lib.divdiffable"></a>

#### divdiffable

```python
def divdiffable(v, u)
```

``v u^{-1}`` if ``u <= v`` with lengths adding (so ``partial_{u}`` can be applied), else ``[]``.

<a id="schubmult.utils.schub_lib.kdown_perms"></a>

#### kdown\_perms

```python
def kdown_perms(perm: Permutation, monoperm: Permutation, p: int,
                k: int) -> list[tuple[Permutation, int, int]]
```

One level of the v-path recursion: all ``(new_perm, pp, sign)`` obtained from ``perm`` by up to ``p``
Bruhat descents through position ``k`` (swapping position ``k-1`` with an untouched position on
either side) such that ``new_perm * monoperm`` has the expected length.

<a id="schubmult.utils.schub_lib.rc_graph_set"></a>

#### rc\_graph\_set

```python
def rc_graph_set(perm)
```

All RC graphs of ``perm`` as ``(row_labels, reduced_word)`` pairs, built recursively with `pull_out_var`.

<a id="schubmult.utils.schub_lib.compute_vpathdicts_cached"></a>

#### compute\_vpathdicts\_cached

```python
@cache
def compute_vpathdicts_cached(th, vmu)
```

Build the v-path dictionaries for theta code ``th`` and target ``vmu``.

Working from the top level down, each level ``i`` maps a permutation to the `kdown_perms`
moves available at that level; the result is re-indexed as
``vpathdicts[i][source] = {(target, degree, sign), ...}`` for forward accumulation.

<a id="schubmult.utils.schub_lib.compute_vpathdicts"></a>

#### compute\_vpathdicts

```python
def compute_vpathdicts(th, vmu)
```

Cached `compute_vpathdicts_cached` accepting a list ``th``.

<a id="schubmult.utils.schub_lib.check_blocks"></a>

#### check\_blocks

```python
def check_blocks(qv, parabolic_index)
```

Parabolic admissibility of a ``q``-exponent vector: over every interval within each block of
consecutive parabolic indices, the sum of `omega` values must be 0 or -1.

<a id="schubmult.utils.schub_lib.reduce_q_coeff"></a>

#### reduce\_q\_coeff

```python
def reduce_q_coeff(u, v, w, qv)
```

One quantum reduction step: find a position where swapping adjacent entries of ``v`` (or ``u``)
together with ``w`` -- adjusting the ``q``-vector when ``w`` had no descent there -- leaves the
quantum LR coefficient unchanged. Returns ``(u, v, w, qv, changed)``.

<a id="schubmult.utils.schub_lib.reduce_q_coeff_u_only"></a>

#### reduce\_q\_coeff\_u\_only

```python
def reduce_q_coeff_u_only(u, v, w, qv)
```

`reduce_q_coeff` restricted to swaps on ``u``.

<a id="schubmult.utils.schub_lib.elem_sym_perms_q"></a>

#### elem\_sym\_perms\_q

```python
def elem_sym_perms_q(orig_perm, p, k, q_var=q_var)
```

Quantum Pieri moves for ``E_p^q(x_1..x_k)``: all ``(perm, degree, q_monomial)`` reachable from
``orig_perm`` by up to ``p`` transpositions ``(i, j)`` with ``i < k <= j`` on still-untouched positions
``i``; a Bruhat *descent* (cyclic quantum move) contributes ``q_{i+1} ... q_j``.

<a id="schubmult.utils.schub_lib.elem_sym_perms_q_op"></a>

#### elem\_sym\_perms\_q\_op

```python
def elem_sym_perms_q_op(orig_perm, p, k, n, q_var=q_var)
```

Adjoint (downward) version of `elem_sym_perms_q` on permutations padded to length ``n``.

<a id="schubmult.utils.schub_lib.elem_sym_perms"></a>

#### elem\_sym\_perms

```python
def elem_sym_perms(orig_perm, p, k)
```

Pieri moves for ``e_p(x_1..x_k)``: all ``(perm, degree)`` reachable from ``orig_perm`` by up to ``p``
Bruhat ascents ``(i, j)`` with ``i < k <= j``, each position ``i`` used at most once (values
swapped in strictly decreasing order to avoid double counting).

<a id="schubmult.utils.schub_lib.elem_sym_perms_groth"></a>

#### elem\_sym\_perms\_groth

```python
def elem_sym_perms_groth(orig_perm, p, k)
```

Grothendieck variant of `elem_sym_perms`: positions in ``i < k`` may be reused, with ``j``
weakly decreasing along the chain.

<a id="schubmult.utils.schub_lib.elem_sym_chains_groth"></a>

#### elem\_sym\_chains\_groth

```python
def elem_sym_chains_groth(orig_perm, p, k)
```

Marked Bruhat chains for the Grothendieck Pieri rule: chains of ascents through position ``k``
together with a marking vector (1 force mark, -1 force unmark, 0 free) enforcing the ordering
conditions on consecutive covers.

<a id="schubmult.utils.schub_lib.elem_sym_positional_perms"></a>

#### elem\_sym\_positional\_perms

```python
def elem_sym_positional_perms(orig_perm, p, *k)
```

Pieri moves for ``e_p`` in an arbitrary set of variable positions ``k`` (1-indexed): Bruhat ascents
swapping an untouched position in ``k`` with a position outside ``k``, returned as
``(perm, degree, sign)`` with sign ``-1`` when the outside position is to the left.

<a id="schubmult.utils.schub_lib.elem_sym_positional_perms_q"></a>

#### elem\_sym\_positional\_perms\_q

```python
def elem_sym_positional_perms_q(orig_perm, p, *k, q_var=q_var)
```

Quantum version of `elem_sym_positional_perms`: returns ``(perm, degree, sign, q_monomial)``.

<a id="schubmult.utils.schub_lib.complete_sym_perms_op"></a>

#### complete\_sym\_perms\_op

```python
def complete_sym_perms_op(orig_perm, p, k)
```

Downward Pieri moves for ``h_p(x_1..x_k)``: ``{perm: [(i, j), ...]}`` recording the chain of Bruhat
descents (position ``i < k`` with an untouched position ``j``).

<a id="schubmult.utils.schub_lib.complete_sym_perms"></a>

#### complete\_sym\_perms

```python
def complete_sym_perms(orig_perm, p, k)
```

Pieri moves for ``h_p(x_1..x_k)``: ``{perm: [(i, j), ...]}`` recording the chain of Bruhat ascents
(position ``i < k`` may repeat; ``j`` untouched).

<a id="schubmult.utils.schub_lib.complete_sym_positional_perms"></a>

#### complete\_sym\_positional\_perms

```python
def complete_sym_positional_perms(orig_perm, p, *k)
```

`complete_sym_perms` for an arbitrary set of positions ``k`` (1-indexed), with signs as in
`elem_sym_positional_perms`; returns ``(perm, degree, sign)`` triples.

<a id="schubmult.utils.schub_lib.elem_sym_perms_op"></a>

#### elem\_sym\_perms\_op

```python
def elem_sym_perms_op(orig_perm, p, k)
```

Downward (Bruhat descent) version of `elem_sym_perms`.

<a id="schubmult.utils.schub_lib.is_split_two"></a>

#### is\_split\_two

```python
def is_split_two(u, v, w)
```

For ``inv(w) - inv(u) == 2``: whether ``v^{-1} w`` is a product of two disjoint cycles; returns
``(True, cycles)`` or ``(False, [])``.

<a id="schubmult.utils.schub_lib.is_coeff_irreducible"></a>

#### is\_coeff\_irreducible

```python
def is_coeff_irreducible(u, v, w)
```

Whether none of the fast paths / reductions applies to the LR triple ``(u, v, w)``.

<a id="schubmult.utils.schub_lib.is_hook"></a>

#### is\_hook

```python
def is_hook(cd)
```

Whether the code ``cd`` is a hook shape: a run of 1's optionally followed by a single larger entry, then zeros.

<a id="schubmult.utils.schub_lib.all_grassmannian_rc_graphs"></a>

#### all\_grassmannian\_rc\_graphs

```python
def all_grassmannian_rc_graphs(n: int, max_inv: int)
```

All RC graphs for Grassmannian permutations generated from partitions.

<a id="schubmult.utils.schub_lib.grassmannian_perms_from_partitions"></a>

#### grassmannian\_perms\_from\_partitions

```python
def grassmannian_perms_from_partitions(n: int,
                                       max_inv: int) -> list[Permutation]
```

Build Grassmannian permutations from partitions: pad to length n, reverse, then uncode.

<a id="schubmult.utils.schub_lib.partitions_with_sum_at_most"></a>

#### partitions\_with\_sum\_at\_most

```python
def partitions_with_sum_at_most(max_sum: int,
                                max_parts: int) -> list[tuple[int, ...]]
```

All partitions (weakly decreasing tuples) with at most max_parts parts and total <= max_sum.

<a id="schubmult.utils.schub_lib.hw_elementary_tensors"></a>

#### hw\_elementary\_tensors

```python
def hw_elementary_tensors(c, bounds=None)
```

Enumerate all sequences (A_1, ..., A_n) with A_i subset of [a_i] and |A_i| = c[i-1],
such that for every j the column sum d_j = #{i : j in A_i} is weakly decreasing in j.

Parameters
----------
c : sequence of nonnegative ints of length n
    A weak composition. Each c[i-1] must satisfy c[i-1] <= a_i.
bounds : sequence of positive ints of length n, optional
    A weakly increasing sequence a_1 <= a_2 <= ... <= a_n of positive integers.
    Defaults to (1, 2, ..., n).

Yields
------
tuple of tuples
    One tuple (A_1, ..., A_n) per valid family; each A_i is a sorted tuple of
    ints in [1, a_i] of length c[i-1].

<a id="schubmult.utils.schub_lib.fff"></a>

#### fff

```python
def fff(chain)
```

Number of forced marks (``1``) in a marked chain from `elem_sym_chains_groth`.

<a id="schubmult.utils.schub_lib.ppp"></a>

#### ppp

```python
def ppp(chain)
```

Number of forced unmarks (``-1``) in a marked chain from `elem_sym_chains_groth`.

<a id="schubmult.utils.schub_lib.groth_pieri_mul"></a>

#### groth\_pieri\_mul

```python
def groth_pieri_mul(perm_dict, p, kk, beta)
```

Grothendieck Pieri rule: multiply the Grothendieck expansion ``perm_dict`` by ``G`` of the
Grassmannian permutation of ``e_p(x_1..x_kk)``. Sums over marked chains from
`elem_sym_chains_groth` with weight ``beta^(length - p) * C(n, k)`` where the free marks are
chosen binomially.

<a id="schubmult.utils.test_utils"></a>

# schubmult.utils.test\_utils

Helpers for the test suite: locating JSON test data and inspecting SymPy/SymEngine expression trees.

<a id="schubmult.utils.test_utils.generate_all"></a>

#### generate\_all

```python
def generate_all(module, filename)
```

Print an import block and ``__all__`` list for the public names defined in ``filename`` (dev helper).

<a id="schubmult.utils.test_utils.get_json"></a>

#### get\_json

```python
def get_json(file: str)
```

Load ``<file>.json`` from the test data directory.

<a id="schubmult.utils.test_utils.load_json_test_names"></a>

#### load\_json\_test\_names

```python
def load_json_test_names(this_dir)
```

Names (without ``.json``) of all test case files in the data subdirectory ``this_dir``.

<a id="schubmult.utils.test_utils.print_args"></a>

#### print\_args

```python
def print_args(poly)
```

Nested string of the argument types of an expression tree (for debugging printing issues).

<a id="schubmult.utils.test_utils.sympify_args"></a>

#### sympify\_args

```python
def sympify_args(poly)
```

Convert a SymPy expression to SymEngine, recursing into ``Mul``/``Pow``/``Add`` when direct
conversion fails.

<a id="schubmult.utils.tuple_utils"></a>

# schubmult.utils.tuple\_utils

Tuple helpers.

<a id="schubmult.utils.tuple_utils.pad_tuple"></a>

#### pad\_tuple

```python
def pad_tuple(tup, length)
```

Pad a tuple-like with trailing zeros up to ``length``.

<a id="schubmult.visualization"></a>

# schubmult.visualization

Visualization tools for Schubert calculus objects.

This module provides functions to visualize combinatorial structures like
pipe dreams from RC graphs.

<a id="schubmult.visualization.draw_pipe_dream"></a>

#### draw\_pipe\_dream

```python
def draw_pipe_dream(rc,
                    max_size=None,
                    title=None,
                    ax=None,
                    flip_horizontal=True,
                    top_labeled=False,
                    show_refs=True)
```

Draw a pipe dream visualization of an RC graph.

In a pipe dream:
- Positions with elements: strands CROSS (go straight through)
- Empty positions: strands AVOID (make 90-degree elbow turns)

Label positioning depends on flip_horizontal and top_labeled:

- flip_horizontal=True, top_labeled=False (default):
Top shows column numbers, right shows permutation output
- flip_horizontal=True, top_labeled=True:
Top shows permutation output, right shows row numbers (1,2,3,...)
- flip_horizontal=False, top_labeled=False:
Left shows permutation output, top shows column numbers
- flip_horizontal=False, top_labeled=True:
Left shows row numbers (1,2,3,...), top shows permutation output

**Arguments**:

- `rc` - RCGraph object to visualize
- `max_size` - Maximum grid size to display (default: determined from permutation)
- `title` - Optional title for the plot
- `ax` - Optional matplotlib axes to draw on (creates new figure if None)
- `flip_horizontal` - If True, reflect horizontally (default: True)
- `top_labeled` - If True, swap which side shows input vs output labels (default: False)
- `show_refs` - If True, display reflection numbers at crossings in green (default: False)
  

**Returns**:

- `tuple` - (fig, ax) matplotlib figure and axes objects
  

**Examples**:

  >>> from schubmult import Permutation, RCGraph
  >>> from schubmult.visualization import draw_pipe_dream
  >>> perm = Permutation([2, 1, 3])
  >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
  >>> fig, ax = draw_pipe_dream(rc, title=f"Pipe Dream for {perm}")
  >>> plt.show()

<a id="schubmult.visualization.draw_pipe_dream_tikz"></a>

#### draw\_pipe\_dream\_tikz

```python
def draw_pipe_dream_tikz(rc,
                         max_size=None,
                         flip_horizontal=True,
                         top_labeled=False,
                         show_refs=False,
                         scale=1.0,
                         outline_rows=None,
                         clip_at_outline=True)
```

Generate TikZ code for a pipe dream visualization of an RC graph.

**Arguments**:

- `rc` - RCGraph object to visualize
- `max_size` - Maximum grid size to display (default: determined from permutation)
- `flip_horizontal` - If True, reflect horizontally (default: True)
- `top_labeled` - If True, swap which side shows input vs output labels (default: False)
- `show_refs` - If True, display reflection numbers at crossings (default: True)
- `scale` - Scale factor for the TikZ picture (default: 1.0)
- `outline_rows` - If provided, draw a thick black outline around this many rows from bottom (default: None)
- `clip_at_outline` - If True and outline_rows is set, clip strands at outline_rows + 1 (default: False)
  

**Returns**:

- `str` - TikZ code as a string
  

**Examples**:

  >>> from schubmult import Permutation, RCGraph
  >>> from schubmult.visualization import draw_pipe_dream_tikz
  >>> perm = Permutation([2, 1, 3])
  >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
  >>> tikz_code = draw_pipe_dream_tikz(rc)
  >>> print(tikz_code)

