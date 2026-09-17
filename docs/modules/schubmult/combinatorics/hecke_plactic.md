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

