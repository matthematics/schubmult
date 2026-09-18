<a id="schubmult.utils.perm_utils"></a>

# schubmult.utils.perm\_utils

Low-level helpers on permutations given as lists/tuples, plus composition and reduced-word utilities.

Most functions here are used by the multiplication kernels in `schubmult.mult` and by
`schubmult.combinatorics.permutation.Permutation`; `add_perm_dict` is the standard way to merge
``{key: coeff}`` expansions.

<a id="schubmult.utils.perm_utils.permtrim_list"></a>

#### permtrim\_list

```python
def permtrim_list(perm)
```

Strip trailing fixed points ``perm[L-1] == L`` from a list in place and return it.

<a id="schubmult.utils.perm_utils.has_bruhat_descent"></a>

#### has\_bruhat\_descent

```python
def has_bruhat_descent(perm, i, j)
```

Check if perm has a Bruhat descent from position i to j.

Optimized version assuming perm is a Permutation object with direct indexing.

<a id="schubmult.utils.perm_utils.count_bruhat"></a>

#### count\_bruhat

```python
def count_bruhat(perm, i, j)
```

Signed length change ``inv(perm * t_{ij}) - inv(perm)`` for the transposition of positions ``i < j``.

<a id="schubmult.utils.perm_utils.has_bruhat_ascent"></a>

#### has\_bruhat\_ascent

```python
def has_bruhat_ascent(perm, i, j)
```

Check if perm has a Bruhat ascent from position i to j.

Optimized version assuming perm is a Permutation object with direct indexing.

<a id="schubmult.utils.perm_utils.omega"></a>

#### omega

```python
def omega(i, qv)
```

``i``-th entry (1-indexed) of the Cartan-matrix image of the q-exponent vector ``qv``:
``2 qv[i] - qv[i-1] - qv[i+1]`` with boundary conventions. Used to convert quantum ``q``
monomials to weights in the parabolic quantum product.

<a id="schubmult.utils.perm_utils.sg"></a>

#### sg

```python
def sg(i, w)
```

1 if ``w`` has a descent at 0-indexed position ``i``, else 0.

<a id="schubmult.utils.perm_utils.count_less_than"></a>

#### count\_less\_than

```python
def count_less_than(arr, val)
```

Number of leading entries of the sorted list ``arr`` that are ``< val``.

<a id="schubmult.utils.perm_utils.artin_sequences"></a>

#### artin\_sequences

```python
def artin_sequences(n)
```

All tuples ``(a_1, ..., a_n)`` with ``0 <= a_i <= n + 1 - i`` (Lehmer codes of ``S_{n+1}``).

<a id="schubmult.utils.perm_utils.weak_compositions"></a>

#### weak\_compositions

```python
def weak_compositions(length, max_degree)
```

All tuples of the given ``length`` with entries in ``0..max_degree``.

<a id="schubmult.utils.perm_utils.is_parabolic"></a>

#### is\_parabolic

```python
def is_parabolic(w, parabolic_index)
```

Whether ``w`` has no descent at any of the (1-indexed) positions in ``parabolic_index``.

<a id="schubmult.utils.perm_utils.add_perm_dict"></a>

#### add\_perm\_dict

```python
def add_perm_dict(d1, d2)
```

Return ``d1 + d2`` as coefficient dicts (keys merged, values added).

<a id="schubmult.utils.perm_utils.add_perm_dict_with_coeff"></a>

#### add\_perm\_dict\_with\_coeff

```python
def add_perm_dict_with_coeff(d1, d2, coeff)
```

Return ``d1 + coeff * d2`` as coefficient dicts.

<a id="schubmult.utils.perm_utils.p_trans"></a>

#### p\_trans

```python
def p_trans(part)
```

Conjugate (transpose) of a partition given as a weakly decreasing list; ``[0]`` for the empty partition.

<a id="schubmult.utils.perm_utils.mu_A"></a>

#### mu\_A

```python
def mu_A(mu, A)
```

The partition whose conjugate consists of the columns of ``mu`` indexed by ``A`` (0-indexed).

<a id="schubmult.utils.perm_utils.get_cycles"></a>

#### get\_cycles

```python
def get_cycles(perm)
```

``perm.get_cycles()``.

<a id="schubmult.utils.perm_utils.old_code"></a>

#### old\_code

```python
def old_code(perm)
```

Lehmer code of a permutation list computed by successive deletion from ``[1..L]``.

<a id="schubmult.utils.perm_utils.cyclic_sort"></a>

#### cyclic\_sort

```python
def cyclic_sort(L)
```

Rotate the list so its maximum is last.

<a id="schubmult.utils.perm_utils.cyclic_sort_min"></a>

#### cyclic\_sort\_min

```python
def cyclic_sort_min(L)
```

Rotate the list so its minimum is first.

<a id="schubmult.utils.perm_utils.h_vector"></a>

#### h\_vector

```python
def h_vector(q_vector)
```

Positions (1-indexed) where the vector strictly increases, up to its first decrease.

<a id="schubmult.utils.perm_utils.l_vector"></a>

#### l\_vector

```python
def l_vector(q_vector)
```

Find l_j = last position where d equals j (where d decreases from j to j-1).

<a id="schubmult.utils.perm_utils.tau_d"></a>

#### tau\_d

```python
def tau_d(d)
```

Partial permutation built from `h_vector`/`l_vector` of ``d`` (``tau[l_i - i] = h_i``), completed by
``Permutation.from_partial``.

<a id="schubmult.utils.perm_utils.phi_d"></a>

#### phi\_d

```python
def phi_d(d)
```

Companion of `tau_d` with the shifted placement ``phi[l_i - 1 - i] = h_i``.

<a id="schubmult.utils.perm_utils.conjugate_weak_composition"></a>

#### conjugate\_weak\_composition

```python
def conjugate_weak_composition(comp)
```

Compute the conjugate of a weak composition.

The conjugate of a weak composition α = (α₁, α₂, ..., αₙ) is the weak composition
β where βⱼ = |{i : αᵢ ≥ j}|, i.e., βⱼ counts how many parts of α are at least j.

This is equivalent to transposing the Ferrers diagram of the composition.

**Arguments**:

- `comp` - A sequence (list, tuple) of non-negative integers representing a weak composition.
  

**Returns**:

  A tuple representing the conjugate weak composition.
  

**Examples**:

  >>> conjugate_weak_composition([3, 1, 0, 2])
  (3, 2, 1)
  >>> conjugate_weak_composition([4, 2, 1])
  (3, 2, 1, 1)
  >>> conjugate_weak_composition([])
  ()
  >>> conjugate_weak_composition([0, 0, 0])
  ()

<a id="schubmult.utils.perm_utils.find_reduced_fail"></a>

#### find\_reduced\_fail

```python
def find_reduced_fail(word, inserted)
```

After changing letter ``inserted`` of a word, find the other position carrying the same root
(the letter whose deletion would make the word reduced again), or ``None``.

<a id="schubmult.utils.perm_utils.is_reduced"></a>

#### is\_reduced

```python
def is_reduced(word)
```

Whether the word of simple reflections is reduced (``inv`` of its product equals its length).

<a id="schubmult.utils.perm_utils.little_bump_pos"></a>

#### little\_bump\_pos

```python
def little_bump_pos(word, index)
```

Little bump at position ``index``: decrement that letter (increment if it is 1), and while the
word is not reduced, repeat at the letter found by `find_reduced_fail`.

<a id="schubmult.utils.perm_utils.little_bump"></a>

#### little\_bump

```python
def little_bump(word, i, j)
```

Little bump of a reduced word at the letter whose right root is the inversion ``(i, j)``.

<a id="schubmult.utils.perm_utils.little_zero"></a>

#### little\_zero

```python
def little_zero(word, length)
```

Repeatedly Little-bump at the last descent until the product's code has fewer than ``length``
entries (Little's map toward a smaller permutation).

