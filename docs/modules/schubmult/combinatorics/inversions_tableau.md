<a id="schubmult.combinatorics.inversions_tableau"></a>

# schubmult.combinatorics.inversions\_tableau

<a id="schubmult.combinatorics.inversions_tableau.InversionsTableau"></a>

## InversionsTableau Objects

```python
class InversionsTableau()
```

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

