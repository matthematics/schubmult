<a id="schubmult.combinatorics.permutation"></a>

# schubmult.combinatorics.permutation

<a id="schubmult.combinatorics.permutation.Permutation"></a>

## Permutation Objects

```python
class Permutation(Printable)
```

Permutation class representing permutations of positive integers.

<a id="schubmult.combinatorics.permutation.Permutation.__truediv__"></a>

#### \_\_truediv\_\_

```python
def __truediv__(other)
```

Returns a tuple of (self, other) if other is a permutation. Intended for skew elements.

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

<a id="schubmult.combinatorics.permutation.Permutation.cycle"></a>

#### cycle

```python
@staticmethod
def cycle(p, q)
```

Construct the cycle permutation used elsewhere in the code.
Kept as a staticmethod on Permutation for call sites like Permutation.cycle(p,q).

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

<a id="schubmult.combinatorics.permutation.Permutation.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

Demazure product

