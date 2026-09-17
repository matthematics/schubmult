<a id="schubmult.combinatorics.increasing_tableau"></a>

# schubmult.combinatorics.increasing\_tableau

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

