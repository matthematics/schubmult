<a id="schubmult.combinatorics.indexed_forests"></a>

# schubmult.combinatorics.indexed\_forests

<a id="schubmult.combinatorics.indexed_forests.IndexedForest"></a>

## IndexedForest Objects

```python
class IndexedForest()
```

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.terminal_nodes"></a>

#### terminal\_nodes

```python
@property
def terminal_nodes()
```

Terminal nodes for trimming, indexed by qdes positions.

In Nadeau--Spink--Tewari notation, qdes(F) is defined from the forest
code c(F) by
    qdes(F) = { i : c_i > 0 and c_{i+1} = 0 }.
We realize these as nodes of index i when present in the support.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.trim_descent_nodes"></a>

#### trim\_descent\_nodes

```python
@property
def trim_descent_nodes()
```

Nodes corresponding to trim descents (qdes positions).

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.trim_descents"></a>

#### trim\_descents

```python
@property
def trim_descents()
```

Trim descents (left terminal set qdes) of the indexed forest.

This matches the code-level criterion used in Nadeau--Spink--Tewari:
    qdes(F) = { i : c_i > 0 and c_{i+1} = 0 },
with 1-based indexing and c_{k}=0 beyond the code length.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.trim_descent"></a>

#### trim\_descent

```python
def trim_descent(index)
```

Trim the forest directly at descent ``index``.

The index must lie in ``self.trim_descents``. If ``c = self.code`` and
``index = i`` (1-based), this performs the inverse of blossoming at i:
decrement ``c_i`` by 1 and delete position ``i+1`` (which is 0 for a
valid trim descent), then rebuild the indexed forest from the resulting
code.

<a id="schubmult.combinatorics.indexed_forests.IndexedForest.is_left_child"></a>

#### is\_left\_child

```python
def is_left_child(descent)
```

Return whether the leaf labeled ``descent`` is a left child.

The label is interpreted as a support index in the indexed forest.
If the labeled node is not a leaf, this returns ``False``.

<a id="schubmult.combinatorics.indexed_forests.draw_indexed_forest"></a>

#### draw\_indexed\_forest

```python
def draw_indexed_forest(forest,
                        save_path=None,
                        show=True,
                        ax=None,
                        dpi=200,
                        node_size=900,
                        support_color="#1f77b4",
                        edge_color="#333333")
```

Draw an indexed forest with support labels.

Parameters
----------
forest : IndexedForest | iterable[Node]
    Forest to draw. This can be the output of `weak_composition_to_indfor`.
save_path : str | None
    If provided, save the figure (e.g. to a `.png` path).
show : bool
    If True and no external axis is provided, display the figure window.
ax : matplotlib.axes.Axes | None
    Optional existing axis to draw into.
dpi : int
    DPI used when creating/saving a new figure.
node_size : int
    Marker size for node circles.
support_color : str
    Color used for support labels.
edge_color : str
    Color used for edges.

Returns
-------
tuple
    `(fig, ax)` for the rendered drawing.

<a id="schubmult.combinatorics.indexed_forests.weak_composition_to_indfor"></a>

#### weak\_composition\_to\_indfor

```python
def weak_composition_to_indfor(c)
```

Converts a weak composition c+ into an indexed forest.
Returns a list of root nodes for the trees in the forest.

<a id="schubmult.combinatorics.indexed_forests.indexed_forest_from_trimming_word"></a>

#### indexed\_forest\_from\_trimming\_word

```python
def indexed_forest_from_trimming_word(trimming_word)
```

Build an indexed forest directly from a trimming word.

If the trimming word is ``(i_1, ..., i_k)``, we reconstruct the forest code
by iterating the inverse blossoming operation

    c <- (c_1, ..., c_{i-1}, c_i + 1, 0, c_{i+1}, ...)

at each index ``i`` from left to right, then convert the resulting weak
composition code into an ``IndexedForest``.

<a id="schubmult.combinatorics.indexed_forests.minimum_n_for_dual_forest"></a>

#### minimum\_n\_for\_dual\_forest

```python
def minimum_n_for_dual_forest(forest_or_code)
```

Smallest n for which the unsigned ~P_F (eq. 4.1) is nonzero.

With L = internal(F): kappa(v) > rho(v), left child < parent (strict),
right child <= parent (weak); we need labels to fit inside [1..n].

<a id="schubmult.combinatorics.indexed_forests.tilde_forest_polynomial"></a>

#### tilde\_forest\_polynomial

```python
def tilde_forest_polynomial(forest_or_code, n=None)
```

Equation (4.1) of arXiv:2306.10939: the unsigned dual forest polynomial.

    ~P_F = sum_{internal(F)-compatible kappa} x^kappa

Returned as an element of FreeAlgebra(WordBasis), keyed by the
composition (exponent tuple) of length n.

<a id="schubmult.combinatorics.indexed_forests.dual_forest_polynomial"></a>

#### dual\_forest\_polynomial

```python
def dual_forest_polynomial(forest_or_code, n)
```

Signed lower-ideal expansion equal to P_F (Thm. 4.1, arXiv:2306.10939).

    P_F = sum_{lower ideals L of F} (-1)^|L| * sum_{L-compatible kappa} x^kappa

Returned as an element of FreeAlgebra(WordBasis). For comparison against
a directly-computed P_F, use the same n and check equality of the
resulting WordBasis dicts.

<a id="schubmult.combinatorics.indexed_forests.build_balanced_tree"></a>

#### build\_balanced\_tree

```python
def build_balanced_tree(labels)
```

Helper to build a tree where the in-order traversal matches the labels.
This creates the 'canonical labeling' referenced in the paper.

<a id="schubmult.combinatorics.indexed_forests.LBS"></a>

## LBS Objects

```python
class LBS(LabeledForest)
```

<a id="schubmult.combinatorics.indexed_forests.LBS.rootlist"></a>

#### rootlist

```python
@property
def rootlist()
```

The finite part of the rootlist: the root labels together with the bare
letters in the gaps of the support.

The genuine rootlist of Nadeau--Tewari Definition 5.6 also contains every
bare letter below and above the support, so it is infinite in both
directions; use :meth:`rootlist_window` to see those tails.

<a id="schubmult.combinatorics.indexed_forests.LBS.rootlist_window"></a>

#### rootlist\_window

```python
def rootlist_window(lo=None, hi=None)
```

Rootlist of Nadeau--Tewari Definition 5.6, truncated to a window.

Consists of the labels of the roots of the trees -- one for each maximal
interval of the support, not one for each node -- together with the bare
letters ``i[0]`` for which neither ``i`` nor ``i-1`` lies in the support.
The bare part is cofinite, so the genuine rootlist is infinite in both
directions; ``lo`` and ``hi`` bound the range of bare letters reported,
defaulting to one step beyond the support on either side.

<a id="schubmult.combinatorics.indexed_forests.LBS.separators"></a>

#### separators

```python
def separators(a, b, lo=None, hi=None)
```

Rootlist elements lying strictly between the letters ``a`` and ``b``.

<a id="schubmult.combinatorics.indexed_forests.LBS.is_separated"></a>

#### is\_separated

```python
def is_separated(a, b)
```

The criterion of Nadeau--Tewari Proposition 5.8 for ``ab <-> ba``.

<a id="schubmult.combinatorics.indexed_forests.omega_reduced_word_from_labelings"></a>

#### omega\_reduced\_word\_from\_labelings

```python
def omega_reduced_word_from_labelings(P: LBS,
                                      Q: DecLabeling) -> tuple[int, ...]
```

Read the reduced-word order from an omega insertion pair (P, Q).

The decreasing labeling Q orders the surviving insertion events; sorting the
forest nodes by Q-label and reading primary labels from P gives the reduced
word represented by the omega insertion output.

<a id="schubmult.combinatorics.indexed_forests.omega_set_valued_compatibility"></a>

#### omega\_set\_valued\_compatibility

```python
def omega_set_valued_compatibility(word_of_pairs: tuple[letterpair, ...])
```

Check set-valued compatibility of a longer sequence via omega insertion.

  Returns a dict with:
      - ok: whether the sequence is set-valued compatible in WCGraph sense,
- omega_reduced_word: reduced word extracted from (P, Q),
- wc_reduced_word: reduced word from WCGraph reduction,
- set_sequence: reduced compatible set sequence from WCGraph.

<a id="schubmult.combinatorics.indexed_forests.omega_is_set_valued_compatible"></a>

#### omega\_is\_set\_valued\_compatible

```python
def omega_is_set_valued_compatible(
        word_of_pairs: tuple[letterpair, ...]) -> bool
```

Boolean convenience wrapper for omega_set_valued_compatibility.

<a id="schubmult.combinatorics.indexed_forests.omega_setvalued_insertion"></a>

#### omega\_setvalued\_insertion

```python
def omega_setvalued_insertion(word, compatible_sequence)
```

Omega insertion for potentially unreduced words with set-valued labels.

This mirrors the reduction logic used by WCGraph: scan a (word, compatible
sequence) pair left-to-right. If appending a letter would be non-reduced,
do not insert a new node; instead merge the compatible label into the set
attached to the corresponding reduced root. Otherwise perform ordinary
insertion on that letterpair.

Returns
-------
dict
    {
      "P": LBS,
      "Q": DecLabeling,
      "reduced_word": tuple[int, ...],
      "set_sequence": tuple[tuple[int, ...], ...],
      "node_sets": dict[int, tuple[int, ...]],
    }
where ``node_sets`` maps forest node index -> merged compatible-label set.

<a id="schubmult.combinatorics.indexed_forests.omega_setvalued_insertion_from_wcgraph"></a>

#### omega\_setvalued\_insertion\_from\_wcgraph

```python
def omega_setvalued_insertion_from_wcgraph(wc_graph)
```

Convenience wrapper: run set-valued omega insertion directly from WCGraph.

<a id="schubmult.combinatorics.indexed_forests.omega_park"></a>

#### omega\_park

```python
def omega_park(word)
```

Parking procedure Omega.

Cars 1, 2, ... arrive successively; car i prefers spot ``word[i]``.
If the preferred spot is empty, the car parks there. Otherwise the
preferred spot lies in a maximal interval [a, b] of occupied spots.
Let j < i be maximal with ``word[j] in [a, b]``; if ``word[i] >= word[j]``
the car parks at b+1, else at a-1.

Returns the set of occupied spots after all cars have parked.

