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

