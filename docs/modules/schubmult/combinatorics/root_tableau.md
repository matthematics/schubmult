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

