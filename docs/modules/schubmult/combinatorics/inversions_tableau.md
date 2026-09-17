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

