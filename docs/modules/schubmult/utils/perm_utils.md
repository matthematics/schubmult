<a id="schubmult.utils.perm_utils"></a>

# schubmult.utils.perm\_utils

<a id="schubmult.utils.perm_utils.has_bruhat_descent"></a>

#### has\_bruhat\_descent

```python
def has_bruhat_descent(perm, i, j)
```

Check if perm has a Bruhat descent from position i to j.

Optimized version assuming perm is a Permutation object with direct indexing.

<a id="schubmult.utils.perm_utils.has_bruhat_ascent"></a>

#### has\_bruhat\_ascent

```python
def has_bruhat_ascent(perm, i, j)
```

Check if perm has a Bruhat ascent from position i to j.

Optimized version assuming perm is a Permutation object with direct indexing.

<a id="schubmult.utils.perm_utils.l_vector"></a>

#### l\_vector

```python
def l_vector(q_vector)
```

Find l_j = last position where d equals j (where d decreases from j to j-1).

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

<a id="schubmult.utils.perm_utils.little_bump_pos"></a>

#### little\_bump\_pos

```python
def little_bump_pos(word, index)
```

Perform a Little bump on a reduced word at the inversion (i, j).

<a id="schubmult.utils.perm_utils.little_bump"></a>

#### little\_bump

```python
def little_bump(word, i, j)
```

Perform a Little bump on a reduced word at the inversion (i, j).

