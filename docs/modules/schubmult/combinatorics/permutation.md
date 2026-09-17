<a id="schubmult.combinatorics.permutation"></a>

# schubmult.combinatorics.permutation

The core indexing object used throughout ``schubmult``: finite permutations.

A `Permutation` wraps a 1-indexed array (``self[i]`` is the value at 0-indexed
position ``i``, extended by the identity ``self[i] = i + 1`` past its window)
and is immutable/hashable/cached (equal permutations of different window
lengths, e.g. ``[2, 1]`` and ``[2, 1, 3]``, compare and hash equal).
Multiplication ``p * q`` is composition of the underlying functions,
``(p * q)[i] = p[q[i] - 1]``; ``~p`` is the inverse. Most of the combinatorics
library is built from a `Permutation`'s Lehmer code (``code``/``trimcode``),
inversions, and Bruhat/weak order.

<a id="schubmult.combinatorics.permutation.Permutation"></a>

## Permutation Objects

```python
class Permutation(Printable)
```

A finite permutation, stored as a 1-indexed array and extended by the identity.

Construct from array form, e.g. ``Permutation([2, 1, 3])``, or via
``uncode(lehmer_code)``. Supports composition (``*``), inversion (``~``),
Bruhat order (``<=``/``bruhat_leq``), and iteration over the window
(``list(perm)``, ``perm[i]`` 0-indexed).

<a id="schubmult.combinatorics.permutation.Permutation.apply"></a>

#### apply

```python
def apply(arr)
```

Permute the elements of ``arr`` by this permutation: return ``[arr[self[i] - 1] for i in range(len(arr))]``.

<a id="schubmult.combinatorics.permutation.Permutation.is_reducible"></a>

#### is\_reducible

```python
@property
def is_reducible()
```

Whether ``self`` splits as a direct sum of two smaller permutations across some fixed point.

<a id="schubmult.combinatorics.permutation.Permutation.reduce"></a>

#### reduce

```python
def reduce(start_spot=1, strict=False)
```

Split ``self`` at a fixed point into ``(perm1, perm2)`` on disjoint blocks of values, or ``None`` if not reducible.

<a id="schubmult.combinatorics.permutation.Permutation.act_root"></a>

#### act\_root

```python
def act_root(a, b)
```

Image of the root/pair ``(a, b)`` (1-indexed positions) under ``self``: ``(self[a-1], self[b-1])``.

<a id="schubmult.combinatorics.permutation.Permutation.one_dominates"></a>

#### one\_dominates

```python
def one_dominates(other)
```

Whether ``self`` one-step dominates ``other`` (one level of the recursive ``dominates`` test).

<a id="schubmult.combinatorics.permutation.Permutation.dominates"></a>

#### dominates

```python
def dominates(other)
```

Whether ``self`` dominates ``other`` in the sense used by the dual Pieri / positivity rules.

<a id="schubmult.combinatorics.permutation.Permutation.antiperm"></a>

#### antiperm

```python
@property
def antiperm()
```

Conjugate of ``self`` by the longest element ``w0``: ``w0 * self * w0``.

<a id="schubmult.combinatorics.permutation.Permutation.weight_coset_decomp"></a>

#### weight\_coset\_decomp

```python
def weight_coset_decomp(dominant_weight)
```

Coset decomposition of ``self`` with respect to the stabilizer of ``dominant_weight``; see ``coset_decomp``.

<a id="schubmult.combinatorics.permutation.Permutation.min_of_weight_coset"></a>

#### min\_of\_weight\_coset

```python
def min_of_weight_coset(dominant_weight)
```

Minimal-length coset representative of ``self`` in the stabilizer of ``dominant_weight``.

<a id="schubmult.combinatorics.permutation.Permutation.max_of_weight_coset"></a>

#### max\_of\_weight\_coset

```python
def max_of_weight_coset(dominant_weight)
```

Maximal-length coset representative of ``self`` in the stabilizer of ``dominant_weight``.

<a id="schubmult.combinatorics.permutation.Permutation.does_demazure_crystal_tensor_decompose"></a>

#### does\_demazure\_crystal\_tensor\_decompose

```python
@staticmethod
def does_demazure_crystal_tensor_decompose(dominant_weight1, low_perm1,
                                           dominant_weight2, low_perm2)
```

Whether the tensor product of the two Demazure crystals (given by their dominant weights and
lowest-weight sorting permutations) decomposes as their `does_demazure_crystal_tensor_decompose`
criterion predicts: the minimal coset representative of ``low_perm1`` in weight ``dominant_weight1``
lies below the longest element of the left descents of the maximal coset representative of
``low_perm2`` in weight ``dominant_weight2``, in Bruhat order.

<a id="schubmult.combinatorics.permutation.Permutation.coset_decomp"></a>

#### coset\_decomp

```python
def coset_decomp(*descs)
```

Decompose ``self = reduced_perm * w_J`` where ``w_J`` lies in the parabolic subgroup generated
by the simple reflections at 1-indexed positions ``descs`` and ``reduced_perm`` is the minimal-length
coset representative (no descents in ``descs``).

**Returns**:

- `tuple` - ``(reduced_perm, w_J)``.

<a id="schubmult.combinatorics.permutation.Permutation.min_coset_rep"></a>

#### min\_coset\_rep

```python
def min_coset_rep(*descs)
```

Minimal-length coset representative of ``self`` for the parabolic subgroup generated at ``descs``.

<a id="schubmult.combinatorics.permutation.Permutation.max_coset_rep"></a>

#### max\_coset\_rep

```python
def max_coset_rep(*descs)
```

Maximal-length coset representative of ``self`` for the parabolic subgroup generated at ``descs``.

<a id="schubmult.combinatorics.permutation.Permutation.longest_element"></a>

#### longest\_element

```python
@classmethod
def longest_element(cls, *descs)
```

Longest element of the parabolic subgroup generated by the simple reflections at 1-indexed positions ``descs``.

<a id="schubmult.combinatorics.permutation.Permutation.w0"></a>

#### w0

```python
@classmethod
def w0(cls, n)
```

The longest element of the symmetric group on ``n`` letters: ``[n, n-1, ..., 1]``.

<a id="schubmult.combinatorics.permutation.Permutation.all_permutations"></a>

#### all\_permutations

```python
@classmethod
@cache
def all_permutations(cls, n)
```

All permutations of ``{1, ..., n}``, as a list of `Permutation`.

<a id="schubmult.combinatorics.permutation.Permutation.parabolic_reduce"></a>

#### parabolic\_reduce

```python
def parabolic_reduce(*descs)
```

Decompose ``self = reduced_perm * w_J`` where ``w_J`` is generated by the simple reflections
*not* in ``descs`` and ``reduced_perm`` has no descents outside ``descs``.

**Returns**:

- `tuple` - ``(reduced_perm, w_J)``.

<a id="schubmult.combinatorics.permutation.Permutation.__truediv__"></a>

#### \_\_truediv\_\_

```python
def __truediv__(other)
```

Returns a tuple of (self, other) if other is a permutation. Intended for skew elements.

<a id="schubmult.combinatorics.permutation.Permutation.fixers"></a>

#### fixers

```python
@staticmethod
def fixers(dominant_weight)
```

1-indexed positions ``i`` where ``dominant_weight[i-1] == dominant_weight[i]``
(the simple reflections fixing ``dominant_weight``, i.e. generating its stabilizer).

<a id="schubmult.combinatorics.permutation.Permutation.ref_product"></a>

#### ref\_product

```python
@classmethod
def ref_product(cls, *args)
```

Product of the simple reflections ``s_a`` for ``a`` in ``args``, applied left to right.

<a id="schubmult.combinatorics.permutation.Permutation.hecke_ref_product"></a>

#### hecke\_ref\_product

```python
@classmethod
def hecke_ref_product(cls, *args)
```

Like ``ref_product``, but each ``s_a`` is applied only if it is a Bruhat ascent
(the 0-Hecke / Demazure product of the simple reflections).

<a id="schubmult.combinatorics.permutation.Permutation.code_word"></a>

#### code\_word

```python
@property
def code_word()
```

A canonical reduced word for ``self``, read off from ``trimcode``.

<a id="schubmult.combinatorics.permutation.Permutation.inverse_code_word"></a>

#### inverse\_code\_word

```python
@property
def inverse_code_word()
```

A canonical reduced word for ``~self``, read off from ``(~self).trimcode``.

<a id="schubmult.combinatorics.permutation.Permutation.root_swap"></a>

#### root\_swap

```python
def root_swap(root)
```

Multiply ``self`` by the reflection swapping the 1-indexed positions ``root = (a, b)``.

<a id="schubmult.combinatorics.permutation.Permutation.reflection"></a>

#### reflection

```python
@classmethod
def reflection(cls, root)
```

The transposition swapping the 1-indexed positions ``root = (a, b)``.

<a id="schubmult.combinatorics.permutation.Permutation.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(index, word=None)
```

The positive root sent to a negative root by the ``index``-th letter of ``word``
(default ``self.code_word``), read from the right (post-multiplied by the remaining suffix).

<a id="schubmult.combinatorics.permutation.Permutation.left_root_at"></a>

#### left\_root\_at

```python
def left_root_at(index, word=None)
```

Like ``right_root_at``, but the root is transported by the prefix of ``word`` before ``index``.

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

<a id="schubmult.combinatorics.permutation.Permutation.commutation_class_of"></a>

#### commutation\_class\_of

```python
@staticmethod
def commutation_class_of(word)
```

All words obtainable from ``word`` by commuting adjacent far-apart letters (``|a - b| >= 2``).

<a id="schubmult.combinatorics.permutation.Permutation.forest_class_of"></a>

#### forest\_class\_of

```python
@staticmethod
def forest_class_of(word)
```

Words reachable from ``word`` by commutation moves that also preserve the indexed-forest
insertion structure (``omega_insertion``); a refinement of ``commutation_class_of``.

<a id="schubmult.combinatorics.permutation.Permutation.code_index_of_index"></a>

#### code\_index\_of\_index

```python
def code_index_of_index(index)
```

The position in ``trimcode`` whose block of ``code_word`` letters contains position ``index``.

<a id="schubmult.combinatorics.permutation.Permutation.cycle"></a>

#### cycle

```python
@staticmethod
def cycle(p, q)
```

Construct the cycle permutation used elsewhere in the code.
Kept as a staticmethod on Permutation for call sites like Permutation.cycle(p,q).

<a id="schubmult.combinatorics.permutation.Permutation.sorting_perm"></a>

#### sorting\_perm

```python
@classmethod
def sorting_perm(cls, itera, reverse=False)
```

The permutation that sorts ``itera`` into (by default) increasing order.

<a id="schubmult.combinatorics.permutation.Permutation.right_act"></a>

#### right\_act

```python
def right_act(lst)
```

Permute the entries of ``lst`` by this permutation, preserving ``lst``'s type (list stays a list).

<a id="schubmult.combinatorics.permutation.Permutation.bruhat_leq"></a>

#### bruhat\_leq

```python
def bruhat_leq(perm, perm2)
```

Whether ``perm <= perm2`` in Bruhat order (equivalently, ``perm``'s tableau criterion
against ``perm2`` on every prefix of their windows).

<a id="schubmult.combinatorics.permutation.Permutation.from_code"></a>

#### from\_code

```python
@classmethod
def from_code(cls, cd)
```

Alias for ``uncode(cd)``: the permutation with Lehmer code ``cd``.

<a id="schubmult.combinatorics.permutation.Permutation.has_pattern"></a>

#### has\_pattern

```python
def has_pattern(pattern)
```

Whether ``self`` contains ``pattern`` (a plain list, not necessarily reduced) as a pattern:
some subsequence of ``self``'s window order-isomorphic to ``pattern``.

<a id="schubmult.combinatorics.permutation.Permutation.zero_indexed_descents"></a>

#### zero\_indexed\_descents

```python
def zero_indexed_descents()
```

0-indexed descent positions: ``i`` such that ``self[i] > self[i+1]``.

<a id="schubmult.combinatorics.permutation.Permutation.descents"></a>

#### descents

```python
def descents(zero_indexed=True)
```

Descent positions of ``self``, 0-indexed by default or 1-indexed if ``zero_indexed=False``.

<a id="schubmult.combinatorics.permutation.Permutation.get_cycles"></a>

#### get\_cycles

```python
def get_cycles(sort_min=False)
```

Cycle decomposition of ``self`` as a list of tuples; ``sort_min`` rotates each cycle to start
at its minimum element instead of the sympy convention.

<a id="schubmult.combinatorics.permutation.Permutation.from_cycles"></a>

#### from\_cycles

```python
@classmethod
def from_cycles(cls, cycle_iter)
```

Build a permutation from a cycle decomposition (sequence of cycles, each a sequence of 1-indexed values).

<a id="schubmult.combinatorics.permutation.Permutation.code"></a>

#### code

```python
@property
def code()
```

Lehmer code of ``self``: ``code[i]`` counts ``j > i`` with ``self[i] > self[j]``.

<a id="schubmult.combinatorics.permutation.Permutation.graph"></a>

#### graph

```python
@property
def graph()
```

The permutation matrix support as a set of 1-indexed pairs ``{(i, self[i])}``.

<a id="schubmult.combinatorics.permutation.Permutation.reduced_with"></a>

#### reduced\_with

```python
def reduced_with(other)
```

Whether ``self * other`` is length-additive: ``inv(self) + inv(other) == inv(self * other)``.

<a id="schubmult.combinatorics.permutation.Permutation.shape"></a>

#### shape

```python
@property
def shape()
```

The Lehmer code sorted into weakly decreasing order (a partition).

<a id="schubmult.combinatorics.permutation.Permutation.inversion_set"></a>

#### inversion\_set

```python
@cached_property
def inversion_set()
```

Set of inversions ``(i, j)`` (1-indexed, ``i < j``) with ``self[i-1] > self[j-1]``.

<a id="schubmult.combinatorics.permutation.Permutation.weak_order_leq"></a>

#### weak\_order\_leq

```python
def weak_order_leq(other)
```

Whether ``self <= other`` in left weak order (``self``'s inversion set is a subset of ``other``'s).

<a id="schubmult.combinatorics.permutation.Permutation.weak_order_meet"></a>

#### weak\_order\_meet

```python
def weak_order_meet(other)
```

Meet (greatest lower bound) of ``self`` and ``other`` in left weak order.

<a id="schubmult.combinatorics.permutation.Permutation.weak_order_join"></a>

#### weak\_order\_join

```python
def weak_order_join(other)
```

Join (least upper bound) of ``self`` and ``other`` in left weak order, via ``w0``-duality with the meet.

<a id="schubmult.combinatorics.permutation.Permutation.diagram"></a>

#### diagram

```python
@property
def diagram()
```

The Rothe diagram of ``self``: cells ``(i, j)`` with ``self[i-1] > j`` and ``(~self)[j-1] > i``.

<a id="schubmult.combinatorics.permutation.Permutation.rothe_diagram"></a>

#### rothe\_diagram

```python
@property
def rothe_diagram()
```

The graph of ``self`` as a set of 1-indexed pairs (see ``diagram`` for the Rothe diagram cells).

<a id="schubmult.combinatorics.permutation.Permutation.max_descent"></a>

#### max\_descent

```python
@cached_property
def max_descent()
```

Number of entries in ``trimcode`` (one past the last nonzero code entry).

<a id="schubmult.combinatorics.permutation.Permutation.maximal_corner"></a>

#### maximal\_corner

```python
@property
def maximal_corner()
```

The maximal corner ``(maxd, end_spot)`` of ``self``'s diagram, used by ``pivots``/``pivot_transition``.

<a id="schubmult.combinatorics.permutation.Permutation.from_partial"></a>

#### from\_partial

```python
@classmethod
def from_partial(cls, partial_perm)
```

Complete a partial assignment (a list with some ``None`` entries) to a full permutation,
filling the gaps with the missing values in increasing order.

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

<a id="schubmult.combinatorics.permutation.Permutation.pad_code"></a>

#### pad\_code

```python
def pad_code(length)
```

``trimcode`` padded with trailing zeros to ``length``.

<a id="schubmult.combinatorics.permutation.Permutation.trimcode"></a>

#### trimcode

```python
@cached_property
def trimcode()
```

Lehmer code truncated to drop trailing zeros (length equals the last descent position).

<a id="schubmult.combinatorics.permutation.Permutation.mul_dominant"></a>

#### mul\_dominant

```python
def mul_dominant()
```

Left factor of ``self`` in its decomposition against the minimal dominant permutation above it.

<a id="schubmult.combinatorics.permutation.Permutation.strict_mul_dominant"></a>

#### strict\_mul\_dominant

```python
def strict_mul_dominant(size=None)
```

Variant of ``mul_dominant`` built from the strict theta (strictly decreasing dominant code).

<a id="schubmult.combinatorics.permutation.Permutation.shiftup"></a>

#### shiftup

```python
def shiftup(k)
```

``self`` with its Lehmer code shifted right by ``k`` (``k`` leading zero code entries prepended).

<a id="schubmult.combinatorics.permutation.Permutation.inv"></a>

#### inv

```python
@cached_property
def inv()
```

Length of ``self``: the number of inversions, i.e. ``sum(self.code)``.

<a id="schubmult.combinatorics.permutation.Permutation.is_dominant"></a>

#### is\_dominant

```python
@property
def is_dominant()
```

Whether ``self`` equals the minimal dominant permutation above it (its code is already weakly decreasing).

<a id="schubmult.combinatorics.permutation.Permutation.is_strict_dominant"></a>

#### is\_strict\_dominant

```python
@property
def is_strict_dominant()
```

Whether ``self`` is dominant with a strictly decreasing ``trimcode``.

<a id="schubmult.combinatorics.permutation.Permutation.is_vexillary"></a>

#### is\_vexillary

```python
@property
def is_vexillary()
```

Whether ``self`` avoids the pattern ``2143`` (equivalently, its Schubert polynomial is a
single Schur polynomial in the ``trimcode``-shape).

<a id="schubmult.combinatorics.permutation.Permutation.swap"></a>

#### swap

```python
def swap(i, j)
```

Multiply ``self`` by the transposition of 0-indexed positions ``i`` and ``j`` (window extended as needed).

<a id="schubmult.combinatorics.permutation.Permutation.rslice"></a>

#### rslice

```python
def rslice(start, stop)
```

Window values at 0-indexed positions ``start`` (inclusive) to ``stop`` (exclusive), extended by fixed points.

<a id="schubmult.combinatorics.permutation.Permutation.__matmul__"></a>

#### \_\_matmul\_\_

```python
def __matmul__(other)
```

Demazure product

<a id="schubmult.combinatorics.permutation.Permutation.pattern_at"></a>

#### pattern\_at

```python
def pattern_at(*indices)
```

The (inverse sorting) pattern induced by ``self`` on the given 1-indexed ``indices``.

<a id="schubmult.combinatorics.permutation.Permutation.minimal_dominant_above"></a>

#### minimal\_dominant\_above

```python
def minimal_dominant_above()
```

The minimal dominant permutation ``>= self`` in Bruhat order: ``uncode(self.theta())``.

<a id="schubmult.combinatorics.permutation.Permutation.foundational_root"></a>

#### foundational\_root

```python
@property
def foundational_root()
```

The pivot pair ``(k, mx)`` marking ``self``'s last descent block and the last position dropping below it.

<a id="schubmult.combinatorics.permutation.Permutation.theta"></a>

#### theta

```python
def theta()
```

Dominant (weakly decreasing) sequence bounding ``self``'s code, used throughout the v-path
multiplication algorithms (see ``schubmult.mult``).

<a id="schubmult.combinatorics.permutation.Permutation.medium_theta"></a>

#### medium\_theta

```python
def medium_theta()
```

Variant of ``theta`` used by the "fast"/merged-layer multiplication kernels.

<a id="schubmult.combinatorics.permutation.Permutation.strict_theta"></a>

#### strict\_theta

```python
def strict_theta()
```

Variant of ``theta`` with strictly decreasing entries (no repeated nonzero layers).

<a id="schubmult.combinatorics.permutation.Permutation.maximal_sortable_below"></a>

#### maximal\_sortable\_below

```python
def maximal_sortable_below()
```

The maximal "sortable" permutation ``<= self`` (code has no gap of more than 1 between
consecutive entries), used by ``mul_sortable``.

<a id="schubmult.combinatorics.permutation.Permutation.mul_sortable"></a>

#### mul\_sortable

```python
def mul_sortable()
```

Left factor of ``self`` against its ``maximal_sortable_below`` decomposition (dual to ``mul_dominant``).

<a id="schubmult.combinatorics.permutation.uncode"></a>

#### uncode

```python
def uncode(cd)
```

The permutation whose Lehmer code is ``cd`` (a list of nonnegative integers).

<a id="schubmult.combinatorics.permutation.permtrim"></a>

#### permtrim

```python
def permtrim(perm)
```

Normalize ``perm`` (a plain array) into a `Permutation` (trims trailing fixed points).

<a id="schubmult.combinatorics.permutation.phi1"></a>

#### phi1

```python
def phi1(u)
```

Drop the first entry of ``(~u).code`` and re-invert: one step of the ``dominates`` recursion.

<a id="schubmult.combinatorics.permutation.split_perms"></a>

#### split\_perms

```python
def split_perms(perms)
```

Split each permutation in ``perms`` (after the first) into two smaller ones across a
reducible point, whenever a valid split point exists; used to normalize a chain of
dominant permutations into minimal reducible pieces.

<a id="schubmult.combinatorics.permutation.s"></a>

#### s

```python
@cache
def s(i)
```

The simple reflection swapping 1-indexed positions ``i`` and ``i + 1``.

