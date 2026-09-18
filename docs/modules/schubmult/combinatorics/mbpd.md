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

