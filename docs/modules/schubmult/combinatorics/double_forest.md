<a id="schubmult.combinatorics.double_forest"></a>

# schubmult.combinatorics.double\_forest

Double forest polynomial construction via the vine subword model.

This module provides a non-script home for the core computational routine
`double_forest_polynomial`, suitable for use by library code.

<a id="schubmult.combinatorics.double_forest.forest_from_code"></a>

#### forest\_from\_code

```python
def forest_from_code(code)
```

Build indexed forest F with c(F)=code via Thompson monoid factorization.

<a id="schubmult.combinatorics.double_forest.forest_qdes"></a>

#### forest\_qdes

```python
def forest_qdes(code)
```

Return the left terminal set (qdes) of the forest with the given code.

qdes(F) = {i : c_i > 0 and c_{i+1} = 0}, using 1-based indexing.

These are the descents of the forest, analogous to descent set of a
permutation (Nadeau-Spink-Tewari, arXiv:2406.01510, §3.2).

<a id="schubmult.combinatorics.double_forest.forest_code_from_trimming_sequence"></a>

#### forest\_code\_from\_trimming\_sequence

```python
def forest_code_from_trimming_sequence(trimming_seq)
```

Return the forest composition sfc(F) for the unique indexed forest F
whose set of trimming sequences contains ``trimming_seq``.

A trimming sequence (i_1, ..., i_k) for F ∈ IndexedForests is defined
recursively: i_k ∈ qdes(F) and (i_1,...,i_{k-1}) ∈ Trim(F/i_k).
(Nadeau-Spink-Tewari, arXiv:2406.01510, Definition 3.8.)

Recovery uses the inverse operation (blossoming):
    sfc(F · i) = (c_1,...,c_{i-1}, c_i+1, 0, c_{i+1}, c_{i+2},...)
i.e. increment position i and insert a 0 immediately after.
Starting from the empty forest and blossoming left-to-right through
(i_1,...,i_k) reconstructs F.

Parameters
----------
trimming_seq : iterable of int
    Sequence of positive integers (1-based leaf positions).

Returns
-------
tuple of int
    The composition sfc(F); trailing zeros are kept to preserve the
    ambient length implied by the sequence.

Examples
--------
>>> forest_code_from_trimming_sequence([1, 1, 2, 4, 7])
(2, 1, 0, 1, 0, 0, 1, 0)
>>> forest_code_from_trimming_sequence([3])
(0, 0, 1, 0)
>>> forest_code_from_trimming_sequence([])
()

<a id="schubmult.combinatorics.double_forest.sylvester_word"></a>

#### sylvester\_word

```python
def sylvester_word(forest)
```

Return one Sylvester word of `forest` via pre-order traversal.

<a id="schubmult.combinatorics.double_forest.sylvester_forest"></a>

#### sylvester\_forest

```python
def sylvester_forest(code, genset, t)
```

Sum of ``polyvalue(genset, t)`` over all RC graphs of ``uncode(code)`` with the given forest weight;
the single (non-double) specialization of `double_sylvester_forest`.

<a id="schubmult.combinatorics.double_forest.double_sylvester_forest"></a>

#### double\_sylvester\_forest

```python
def double_sylvester_forest(code, genset, t)
```

Double (equivariant) forest polynomial, computed by pairing RC graphs of ``u`` and ``v`` from
the double Schubert expansion ``Sx([]) * DSx(perm, "t")`` whose merged vine diagram matches the
principal RC graph's omega-invariant target.

<a id="schubmult.combinatorics.double_forest.canonical_forest_from_word"></a>

#### canonical\_forest\_from\_word

```python
def canonical_forest_from_word(word)
```

Canonical forest form used to test Sylvester-equivalence of words.

<a id="schubmult.combinatorics.double_forest.long_word"></a>

#### long\_word

```python
def long_word(n)
```

Build omega_tilde_[n] as a list of letter records.

<a id="schubmult.combinatorics.double_forest.letter_weight"></a>

#### letter\_weight

```python
def letter_weight(letter, x_gen, t_gen)
```

Weight map in the vine model. x_gen, t_gen are 1-indexed subscriptables.

<a id="schubmult.combinatorics.double_forest.double_forest_polynomial"></a>

#### double\_forest\_polynomial

```python
def double_forest_polynomial(code, x_gen, t_gen, n=None)
```

Compute P_F(x;t) for indexed forest F with code c(F)=code.

<a id="schubmult.combinatorics.double_forest.reflection_subword_polynomial"></a>

#### reflection\_subword\_polynomial

```python
def reflection_subword_polynomial(code, x_gen, t_gen, n=None, mode="forest")
```

Vine subword model from arXiv:2504.15234 (Bergeron-Gagnon-Nadeau-Spink-Tewari).

Iterates subwords pi of the long word `long_word(n)` of size
    sz = sum(code) = ell(uncode(code))
and sums wt(pi) under one of two filters on the *value sequence*
(treating barred and unbarred letters by their value only):

  mode='forest'   :  values must be Sylvester-equivalent to sylvester_word(F),
                     i.e. lie in Syl(F).  Equals P_F(x;t)  (Theorem 5.1).
  mode='schubert' :  values must be a reduced word for perm=uncode(code),
                     i.e. lie in Red(perm).  Equals S_perm(x;t) (Theorem 6.1).

Weights (Sylvester column convention, paper p.19):
    wt(j^(k))      = x_k - t_j
    wt(barred j^(k)) = t_j - t_k
realised by `letter_weight` with bracket-indexed `x_gen`, `t_gen`.

<a id="schubmult.combinatorics.double_forest.reflection_forest_polynomial"></a>

#### reflection\_forest\_polynomial

```python
def reflection_forest_polynomial(code, x_gen, t_gen, n=None, weight_rule=None)
```

Reflection-style forest model from the FULL long word (barred + unbarred).

Iterates subwords of `long_word(n)` (same enumeration used by
`double_forest_polynomial`). A subword is kept iff:
  (a) its value sequence is a reduced expression for the same
      permutation as `sylvester_word(F)`, AND
  (b) its omega-insertion P-symbol (Nadeau-Tewari arXiv:2306.10939 §5.1)
      equals that of `sylvester_word(F)` (single Omega-class).

Default weight is `letter_weight` (paper convention). Override with
`weight_rule(letter, position_in_subword, picked_letters, x_gen, t_gen,
             sylvester_letters)` -> sympy expression, where
    letter            = picked long_word letter dict
    position_in_subword = 0-indexed position in the picked subword
    picked_letters    = tuple of all picked long_word letter dicts
    sylvester_letters = tuple of long_word letters of the canonical
                        Sylvester subword that gave `target`

<a id="schubmult.combinatorics.double_forest.debug_compare_models"></a>

#### debug\_compare\_models

```python
def debug_compare_models(code, x_gen, t_gen, n=None)
```

Print, for `code`, the subwords that survive in:

(A) the alphabet vine model (double_forest_polynomial), and
(B) the reflection-alphabet model with the omega-insertion filter
    using the principal-RC reversed perm-word as target,

so we can stare at the symmetric difference.

