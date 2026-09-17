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

