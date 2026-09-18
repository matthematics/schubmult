<a id="schubmult.rings.free_algebra.schubert_basis"></a>

# schubmult.rings.free\_algebra.schubert\_basis

`SchubertBasis`: the free-algebra basis dual to Schubert polynomials.

A key is ``(perm, numvars)``: the element dual to ``S_perm`` viewed as a polynomial in
exactly ``numvars`` variables (so ``numvars >= max_descent(perm)``). Under the
word/monomial pairing this is the ``SchubertPolyBasis`` of the polynomial algebra.

The product is the separated-descents product (`SeparatedDescentsRing`): ``(u, p) * (v, q)``
places ``u`` in the first ``p`` variables and ``v`` in the next ``q``, giving a ``(w, p + q)``
expansion. ``transition_word`` expands a key into words via the SEM (elementary symmetric)
factorization of ``S_perm``, and the other ``transition_*`` methods reach the remaining
bases either directly or by way of the word basis. ``ASx`` is the standard instance.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis"></a>

## SchubertBasis Objects

```python
class SchubertBasis(FreeAlgebraBasis)
```

Free-algebra basis dual to Schubert polynomials; keys are ``(Permutation, numvars)``.

See the module docstring. Products go through the separated-descents Schubert ring,
and transitions to the word basis go through elementary symmetric function
decompositions.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Whether ``x`` is ``(perm,)`` or ``(perm, numvars)`` with ``perm`` a permutation/list/tuple.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

The key ``(rc_graph.perm, len(rc_graph))``: an RC graph's permutation with its row count as ``numvars``.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize to ``(Permutation, numvars)``; if ``numvars`` is omitted it defaults to the last descent.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Separated-descents product: ``(u, p) * (v, q)`` with ``u`` in the first ``p`` variables and
``v`` in the next ``q``, computed in `SeparatedDescentsRing`.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

The skew element ``S_w / S_u`` in ``n`` variables: the dual of multiplying by ``S_u``,
computed with the descent-side kernel ``schubmult_py_down`` and truncated to permutations
fitting in ``n`` variables.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Coproduct of ``(perm, numvars)`` (dual to polynomial multiplication): expand to words,
apply the word coproduct, and convert each tensor factor back to Schubert keys.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in the `GrothendieckBasis`, by taking the co-BPD of every RC graph
of ``perm * w0`` and collecting the resulting permutations.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schubert_schur"></a>

#### transition\_schubert\_schur

```python
@classmethod
def transition_schubert_schur(cls, *x)
```

Expand ``(perm, numvars)`` in the `SchubertSchurBasis`: split off the variables beyond
``numvars`` via a Schubert coproduct against a dominant permutation, yielding
``(partition, perm', numvars)`` keys.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schur_elementary"></a>

#### transition\_schur\_elementary

```python
@classmethod
def transition_schur_elementary(cls, *x)
```

Expand ``(perm, numvars)`` in the `SchurElementaryBasis` (a word-like tuple paired with a partition).

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_elementary"></a>

#### transition\_elementary

```python
@classmethod
def transition_elementary(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in the `ElementaryBasis`: read the monomials of ``S_{perm * w0}``
and complement each exponent against the staircase to get elementary-symmetric indices.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_separated_descents"></a>

#### transition\_separated\_descents

```python
@classmethod
def transition_separated_descents(cls, k, *x)
```

Expand ``(perm, numvars)`` in the level-``k`` `SeparatedDescentsBasis` via a Schubert coproduct
splitting the last ``k - 1`` variables, yielding ``(perm_left, perm_right, numvars)`` keys.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, perm, n)
```

Expand ``(perm, n)`` in the `JBasis`: a code with no zeros is already a J key; leading zeros are
peeled off (each contributing a factor ``t``, currently ``1``), and anything else goes via words.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

``SchubertPolyBasis``: Schubert polynomials are the dual basis under the word/monomial pairing.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition"></a>

#### transition

```python
@classmethod
@cache
def transition(cls, other_basis)
```

Return the key -> ``{key: coeff}`` function into ``other_basis``.

Direct routes exist for the word, elementary, Schubert-Schur, Schur-elementary,
composition-Schubert, separated-descents, and Grothendieck bases; everything else
is reached by going through the word basis first.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.old_transition_word"></a>

#### old\_transition\_word

```python
@classmethod
@cache
def old_transition_word(cls, perm, numvars)
```

Transition to WordBasis via SEM basis (legacy implementation).

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, perm, numvars)
```

Expand ``(perm, numvars)`` in the word basis.

Multiplies ``perm`` by the inverse of the dominant permutation with code
``(inv + numvars, ..., inv + 1)`` (``inv = perm.inv``) and writes the result in the SEM
basis with a custom ``elem_func`` that records each factor ``e_p(x_1..x_k)`` as the word
with ``k - p`` in position ``numvars - k + inv``. Products of such words in the monomial
polynomial algebra add letterwise, so the resulting polynomial *is* the word expansion.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Display symbol for ``(perm, numvars)`` (the separated-descents ``SepDescSchubPoly`` form).

