<a id="schubmult.combinatorics.nilplactic"></a>

# schubmult.combinatorics.nilplactic

<a id="schubmult.combinatorics.nilplactic.NilPlactic"></a>

## NilPlactic Objects

```python
class NilPlactic(Plactic)
```

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

