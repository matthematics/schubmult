<a id="schubmult"></a>

# schubmult

<a id="schubmult.__getattr__"></a>

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

Lazy import exported names. This allows `import schubmult` to succeed
even if some optional dependencies are missing.

<a id="schubmult.__dir__"></a>

#### \_\_dir\_\_

```python
def __dir__()
```

Include lazily-exported names in dir()

<a id="schubmult.visualization"></a>

# schubmult.visualization

Visualization tools for Schubert calculus objects.

This module provides functions to visualize combinatorial structures like
pipe dreams from RC graphs.

<a id="schubmult.visualization.draw_pipe_dream"></a>

#### draw\_pipe\_dream

```python
def draw_pipe_dream(rc,
                    max_size=None,
                    title=None,
                    ax=None,
                    flip_horizontal=True,
                    top_labeled=False,
                    show_refs=True)
```

Draw a pipe dream visualization of an RC graph.

In a pipe dream:
- Positions with elements: strands CROSS (go straight through)
- Empty positions: strands AVOID (make 90-degree elbow turns)

Label positioning depends on flip_horizontal and top_labeled:

- flip_horizontal=True, top_labeled=False (default):
Top shows column numbers, right shows permutation output
- flip_horizontal=True, top_labeled=True:
Top shows permutation output, right shows row numbers (1,2,3,...)
- flip_horizontal=False, top_labeled=False:
Left shows permutation output, top shows column numbers
- flip_horizontal=False, top_labeled=True:
Left shows row numbers (1,2,3,...), top shows permutation output

**Arguments**:

- `rc` - RCGraph object to visualize
- `max_size` - Maximum grid size to display (default: determined from permutation)
- `title` - Optional title for the plot
- `ax` - Optional matplotlib axes to draw on (creates new figure if None)
- `flip_horizontal` - If True, reflect horizontally (default: True)
- `top_labeled` - If True, swap which side shows input vs output labels (default: False)
- `show_refs` - If True, display reflection numbers at crossings in green (default: False)
  

**Returns**:

- `tuple` - (fig, ax) matplotlib figure and axes objects
  

**Examples**:

  >>> from schubmult import Permutation, RCGraph
  >>> from schubmult.visualization import draw_pipe_dream
  >>> perm = Permutation([2, 1, 3])
  >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
  >>> fig, ax = draw_pipe_dream(rc, title=f"Pipe Dream for {perm}")
  >>> plt.show()

<a id="schubmult.visualization.draw_pipe_dream_tikz"></a>

#### draw\_pipe\_dream\_tikz

```python
def draw_pipe_dream_tikz(rc,
                         max_size=None,
                         flip_horizontal=True,
                         top_labeled=False,
                         show_refs=False,
                         scale=1.0,
                         outline_rows=None,
                         clip_at_outline=True)
```

Generate TikZ code for a pipe dream visualization of an RC graph.

**Arguments**:

- `rc` - RCGraph object to visualize
- `max_size` - Maximum grid size to display (default: determined from permutation)
- `flip_horizontal` - If True, reflect horizontally (default: True)
- `top_labeled` - If True, swap which side shows input vs output labels (default: False)
- `show_refs` - If True, display reflection numbers at crossings (default: True)
- `scale` - Scale factor for the TikZ picture (default: 1.0)
- `outline_rows` - If provided, draw a thick black outline around this many rows from bottom (default: None)
- `clip_at_outline` - If True and outline_rows is set, clip strands at outline_rows + 1 (default: False)
  

**Returns**:

- `str` - TikZ code as a string
  

**Examples**:

  >>> from schubmult import Permutation, RCGraph
  >>> from schubmult.visualization import draw_pipe_dream_tikz
  >>> perm = Permutation([2, 1, 3])
  >>> rc = list(RCGraph.all_rc_graphs(perm, 2))[0]
  >>> tikz_code = draw_pipe_dream_tikz(rc)
  >>> print(tikz_code)

<a id="schubmult.rings.quasisymmetric_functions"></a>

# schubmult.rings.quasisymmetric\_functions

<a id="schubmult.rings.quasisymmetric_functions.stuffle"></a>

#### stuffle

```python
def stuffle(alpha, beta)
```

Computes the stuffle product of two compositions alpha and beta.
Returns a dictionary where keys are resulting compositions (tuples)
and values are their coefficients.

<a id="schubmult.rings.quasisymmetric_functions.quasi_schur_to_monomial"></a>

#### quasi\_schur\_to\_monomial

```python
def quasi_schur_to_monomial(comp)
```

Computes the quasi-Schur function for composition comp in the monomial basis.
Returns a dictionary where keys are compositions (tuples) and values are coefficients.

Uses the standard composition tableau definition: sum over all descent compositions
of standard composition tableaux of the given shape.

<a id="schubmult.rings.quasisymmetric_functions.QSymElement"></a>

## QSymElement Objects

```python
class QSymElement(BaseSchubertElement)
```

<a id="schubmult.rings.quasisymmetric_functions.QSymElement.expand"></a>

#### expand

```python
def expand(num_vars)
```

Expand the quasi-symmetric function in the given number of variables.

<a id="schubmult.rings.quasisymmetric_functions.QSym"></a>

## QSym Objects

```python
class QSym(BaseSchubertRing)
```

<a id="schubmult.rings.quasisymmetric_functions.QSym.quasi_schur"></a>

#### quasi\_schur

```python
def quasi_schur(*comp)
```

Returns the quasi-Schur function for the given composition
expressed in the monomial basis.

**Arguments**:

- ```*comp``` - A composition (tuple or sequence of positive integers)
  

**Returns**:

  QSymElement representing the quasi-Schur function in monomial basis
  

**Example**:

  >>> QS = QSym()
  >>> QS.quasi_schur(2, 1)  # quasi-Schur function for composition (2,1)

<a id="schubmult.rings.free_algebra.schur_elementary_basis"></a>

# schubmult.rings.free\_algebra.schur\_elementary\_basis

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis"></a>

## SchurElementaryBasis Objects

```python
class SchurElementaryBasis(FreeAlgebraBasis)
```

Schur-Elementary basis of the free algebra.

Keys are ``(tuple, tuple)`` pairs encoding products of
a standard elementary monomial (first tuple) and a Schur polynomial (second tuple, partition
of length precisely len(first_tuple) + 1 in increasing order, with zeros at the beginning if
needed).  For example, the key ``((1, 2), (0, 1, 3))`` corresponds to the product of the elementary monomial

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(list/tuple, list/tuple)`` pair.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, tuple)`` key.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert-Schur key via the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, elem_tup, lambd)
```

Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, elem_tup, lambd)
```

Transition a Schubert-Schur key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from SchurElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schur_elementary_basis.SchurElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``SE``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.free_algebra_basis"></a>

# schubmult.rings.free\_algebra.free\_algebra\_basis

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis"></a>

## FreeAlgebraBasis Objects

```python
class FreeAlgebraBasis()
```

Abstract base class for free algebra bases.

Subclasses define how keys are represented, how products and coproducts
are computed, and how to transition between bases.  Default implementations
delegate through the :class:`WordBasis` via ``compose_transition``.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid key for this basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Convert an RC graph to a basis-keyed dict.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a canonical key for this basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a function mapping keys of this basis to dicts in *other_basis*.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, key)
```

Return the display symbol for *key*.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.compose_transition"></a>

#### compose\_transition

```python
@classmethod
def compose_transition(cls, tkeyfunc, output)
```

Apply a key-level transition function to each key in *output*.

For each ``(key, v)`` in *output*, expands ``tkeyfunc(key)`` and
accumulates the results weighted by *v*.

**Arguments**:

- `tkeyfunc` - A function mapping a key to a ``{key: coeff}`` dict.
- `output` - A ``{key: coeff}`` dict to transform.
  

**Returns**:

  A merged ``{key: coeff}`` dict in the target basis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual basis class (for polynomial algebra pairing).

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.change_tensor_basis"></a>

#### change\_tensor\_basis

```python
@classmethod
def change_tensor_basis(cls, tensor_elem, basis1, basis2)
```

Change the bases of both factors of a tensor element.

**Arguments**:

- `tensor_elem` - An element of a tensor product ring.
- `basis1` - Target basis for the left factor.
- `basis2` - Target basis for the right factor.
  

**Returns**:

  The tensor element re-expressed in the new bases.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of *key* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key)
```

Compute the bar-coproduct of *key* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product of two keys by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Inject *key2* into *key1* at position *i* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.prefix"></a>

#### prefix

```python
@classmethod
def prefix(cls, key, length, coeff=S.One)
```

Extract a prefix of *length* letters by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.suffix"></a>

#### suffix

```python
@classmethod
def suffix(cls, key, length, coeff=S.One)
```

Extract a suffix of *length* letters by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.free_algebra_basis.FreeAlgebraBasis.interval"></a>

#### interval

```python
@classmethod
def interval(cls, key, start, stop, coeff=S.One)
```

Extract a subword from *start* to *stop* by delegating through WordBasis.

<a id="schubmult.rings.free_algebra.separated_descents_basis"></a>

# schubmult.rings.free\_algebra.separated\_descents\_basis

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis"></a>

## \_SeparatedDescentsBasis Objects

```python
class _SeparatedDescentsBasis(FreeAlgebraBasis)
```

Separated descents basis of the free algebra (parameterized by level *k*).

Keys are ``(Permutation, Permutation, int)`` triples representing a
factorization into descents above and below a cutoff level *k*.
Instances are created by the :func:`SeparatedDescentsBasis` factory.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid separated descents key.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(Permutation, Permutation, int)`` key.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two separated descents keys via the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, perm0, perm1, numvars)
```

Transition a separated descents key to the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, perm0, perm1, n)
```

Transition a separated descents key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from this separated descents basis to *other_basis*.

<a id="schubmult.rings.free_algebra.separated_descents_basis._SeparatedDescentsBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``SepDesc<k>``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.separated_descents_basis.SeparatedDescentsBasis"></a>

#### SeparatedDescentsBasis

```python
def SeparatedDescentsBasis(k)
```

Factory that creates a separated descents basis class for level *k*.

<a id="schubmult.rings.free_algebra.composition_schubert_basis"></a>

# schubmult.rings.free\_algebra.composition\_schubert\_basis

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis"></a>

## CompositionSchubertBasis Objects

```python
class CompositionSchubertBasis(FreeAlgebraBasis)
```

Schubert basis indexed by padded trimcode compositions.

A key is a composition ``c`` and corresponds to Schubert key
``(uncode(c), len(c))``.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list (composition).

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.as_schubert_key"></a>

#### as\_schubert\_key

```python
@classmethod
def as_schubert_key(cls, key)
```

Convert a composition key to a Schubert key ``(Permutation, length)``.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, key)
```

Normalize a key to a composition tuple.

Accepts either a Schubert key ``(Permutation, int)`` or a raw tuple.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the composition key for the given RC graph.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Inject *key2* into *key1* at position *i* using Schubert multiplication.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two composition keys by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
def coproduct(cls, key)
```

Compute the coproduct by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
def bcoproduct(cls, key)
```

Compute the bar-coproduct by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

Compute the skew element by delegating to SchubertBasis.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual basis (delegates to SchubertBasis).

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from CompositionSchubertBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.composition_schubert_basis.CompositionSchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``CompSchub``-labelled display object for the composition key *k*.

<a id="schubmult.rings.free_algebra.schubert_basis"></a>

# schubmult.rings.free\_algebra.schubert\_basis

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis"></a>

## SchubertBasis Objects

```python
class SchubertBasis(FreeAlgebraBasis)
```

Schubert basis of the free algebra.

Keys are ``(Permutation, int)`` pairs where the integer records the
number of variables.  Products use the separated-descents Schubert
ring, and transitions to the word basis go through elementary
symmetric function decompositions.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a valid Schubert basis key.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the Schubert key ``(perm, length)`` for the given RC graph.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(Permutation, numvars)`` key.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two Schubert basis keys via the separated-descents ring.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.skew_element"></a>

#### skew\_element

```python
@classmethod
def skew_element(cls, w, u, n)
```

Compute the skew Schubert element S_w / S_u truncated to *n* variables.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert key in the tensor ring.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, perm, numvars)
```

Transition a Schubert key to the Grothendieck basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schubert_schur"></a>

#### transition\_schubert\_schur

```python
@classmethod
def transition_schubert_schur(cls, *x)
```

Transition a Schubert key to the Schubert-Schur basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_schur_elementary"></a>

#### transition\_schur\_elementary

```python
@classmethod
def transition_schur_elementary(cls, *x)
```

Transition a Schubert key to the Schur-Elementary basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_elementary"></a>

#### transition\_elementary

```python
@classmethod
def transition_elementary(cls, perm, numvars)
```

Transition a Schubert key to the elementary basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_separated_descents"></a>

#### transition\_separated\_descents

```python
@classmethod
def transition_separated_descents(cls, k, *x)
```

Transition a Schubert key to the separated descents basis of level *k*.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, perm, n)
```

Transition a Schubert key ``(perm, n)`` to the J basis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the SchubertPolyBasis as the dual of SchubertBasis.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.transition"></a>

#### transition

```python
@classmethod
@cache
def transition(cls, other_basis)
```

Return a transition function from SchubertBasis to *other_basis*.

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

Transition a Schubert key to the word basis via SEM factorization.

<a id="schubmult.rings.free_algebra.schubert_basis.SchubertBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a Xi-notation display object for key *k*.

<a id="schubmult.rings.free_algebra.grothendieck_basis"></a>

# schubmult.rings.free\_algebra.grothendieck\_basis

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis"></a>

## GrothendieckBasis Objects

```python
class GrothendieckBasis(FreeAlgebraBasis)
```

Grothendieck basis for FreeAlgebra.

Keys match the Schubert basis shape: ``(Permutation, numvars)``.
The deformation parameter ``beta`` is a class variable so the class can be
passed directly to ``FreeAlgebra`` without requiring a basis constructor.

<a id="schubmult.rings.free_algebra.grothendieck_basis.GrothendieckBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.schubert_schur_basis"></a>

# schubmult.rings.free\_algebra.schubert\_schur\_basis

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis"></a>

## SchubertSchurBasis Objects

```python
class SchubertSchurBasis(FreeAlgebraBasis)
```

Schubert-Schur basis of the free algebra.

Keys are ``(partition_tuple, Permutation)`` pairs encoding products of
a Grassmannian Schubert polynomial with a Schur polynomial.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(list/tuple, Permutation/list/tuple)`` pair.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, Permutation)`` key.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key)
```

Compute the coproduct of a Schubert-Schur key via the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, lambd, perm)
```

Transition a Schubert-Schur key ``(lambda, perm)`` to the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, lambd, perm)
```

Transition a Schubert-Schur key to the word basis via the Schubert basis.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from SchubertSchurBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.schubert_schur_basis.SchubertSchurBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``SS``-prefixed symbol for key *k*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis"></a>

# schubmult.rings.free\_algebra.monomial\_slide\_basis

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis"></a>

## MonomialSlideBasis Objects

```python
class MonomialSlideBasis(FreeAlgebraBasis)
```

Monomial slide basis of the free algebra.

Keys are weak composition tuples. Transitions use monomial slide
polynomial expansions and coarsenings of compositions.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``MS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
@classmethod
@cache
def transition_fundamental_slide(cls, key)
```

Transition a monomial slide key to the fundamental slide basis.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from MonomialSlideBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.monomial_slide_basis.MonomialSlideBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, key)
```

Transition a monomial slide key to the word basis.

<a id="schubmult.rings.free_algebra.word_basis"></a>

# schubmult.rings.free\_algebra.word\_basis

<a id="schubmult.rings.free_algebra.word_basis.WordBasis"></a>

## WordBasis Objects

```python
class WordBasis(FreeAlgebraBasis)
```

Word basis of the free algebra.

Keys are tuples of nonnegative integers representing words. This is the
fundamental basis through which all other bases perform their operations
via basis transitions.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Return the length vector of the RC graph as a word key.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Concatenate two words.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.inject"></a>

#### inject

```python
@classmethod
def inject(cls, key1, i, key2, coeff=S.One)
```

Insert *key2* into *key1* at position *i*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.prefix"></a>

#### prefix

```python
@classmethod
def prefix(cls, key, length, coeff=S.One)
```

Return the first *length* letters of *key*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.suffix"></a>

#### suffix

```python
@classmethod
def suffix(cls, key, length, coeff=S.One)
```

Return the last *length* letters of *key*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.interval"></a>

#### interval

```python
@classmethod
def interval(cls, key, start, stop, coeff=S.One)
```

Return the subword ``key[start:stop]``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.coproduct"></a>

#### coproduct

```python
@classmethod
@cache
def coproduct(cls, key, coeff=S.One)
```

Compute the additive coproduct of a word.

Decomposes each letter into all (i, key-i) splittings and combines
via a divide-and-conquer tensor product.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.bcoproduct"></a>

#### bcoproduct

```python
@classmethod
@cache
def bcoproduct(cls, key, coeff=S.One)
```

Compute the bar-coproduct of a word.

Like :meth:`coproduct` but drops empty factors (zeros map to
the empty tuple).

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.try_internal_product"></a>

#### try\_internal\_product

```python
@classmethod
def try_internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product via integer matrices (requires SageMath).

Uses shifted keys (incremented by 1) with ``IntegerMatrices``.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.internal_product"></a>

#### internal\_product

```python
@classmethod
def internal_product(cls, key1, key2, coeff=S.One)
```

Compute the internal product of two words via integer matrices (requires SageMath).

Words must not contain zeros. Returns the dict of result
words weighted by *coeff*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a bracket-notation symbol like ``[210]`` for the word *k*.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Expand a word tuple into the single Schubert basis via divide-and-conquer.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.jbasis_tup_expand"></a>

#### jbasis\_tup\_expand

```python
@staticmethod
@cache
def jbasis_tup_expand(tup)
```

Expand a word tuple into the Z basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a word key to the Schubert basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_jbasis"></a>

#### transition\_jbasis

```python
@classmethod
def transition_jbasis(cls, key)
```

Transition a word key to the J basis via Pieri products.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_jtbasis"></a>

#### transition\_jtbasis

```python
@classmethod
def transition_jtbasis(cls, key)
```

Transition a word key to the JT basis via normalization.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_forest"></a>

#### transition\_forest

```python
@classmethod
def transition_forest(cls, key)
```

Transition a word key to the forest basis via RC graph enumeration.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_grove"></a>

#### transition\_grove

```python
@classmethod
def transition_grove(cls, key)
```

Transition a word key to the grove basis via WC graph enumeration.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the MonomialBasis as the dual of WordBasis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_monomial_slide"></a>

#### transition\_monomial\_slide

```python
@classmethod
def transition_monomial_slide(cls, key)
```

Transition a word key to the monomial slide basis.

<a id="schubmult.rings.free_algebra.word_basis.WordBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
@cache
def transition_grothendieck(cls, key)
```

Transition a word key (composition) to the Grothendieck basis.

Coefficient of ``G_w`` is the number of WC graphs of permutation ``w``
and weight ``key``, multiplied by ``beta^(|key|-inv(w))``.

<a id="schubmult.rings.free_algebra.elementary_basis"></a>

# schubmult.rings.free\_algebra.elementary\_basis

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis"></a>

## ElementaryBasis Objects

```python
class ElementaryBasis(FreeAlgebraBasis)
```

Elementary symmetric function basis of the free algebra.

Keys are ``(tuple, int)`` pairs where the tuple encodes an elementary
symmetric function composition and the integer is the number of variables.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a ``(tuple/list, int)`` pair.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* into a ``(tuple, int)`` key.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, tup, numvars)
```

Transition an elementary key to the Schubert basis.

<a id="schubmult.rings.free_algebra.elementary_basis.ElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``Elem``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.forest_basis"></a>

# schubmult.rings.free\_algebra.forest\_basis

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis"></a>

## ForestBasis Objects

```python
class ForestBasis(FreeAlgebraBasis)
```

Forest basis of the free algebra.

Keys are tuples representing indexed-forest weight vectors.
Transitions to the Schubert basis use RC graph enumeration
filtered by forest weight.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``ForestDual``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the ForestPolyBasis as the dual of ForestBasis.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a forest key to the Schubert basis via RC graph enumeration.

<a id="schubmult.rings.free_algebra.forest_basis.ForestBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ForestBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.grove_basis"></a>

# schubmult.rings.free\_algebra.grove\_basis

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis"></a>

## GroveBasis Objects

```python
class GroveBasis(FreeAlgebraBasis)
```

Grove basis of the free algebra.

Keys are tuples representing indexed-grove weight vectors.
Transitions to the Grothendieck basis use RC graph enumeration
filtered by grove weight.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``GroveDual``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the GrovePolyBasis as the dual of GroveBasis.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a grove key to the Grothendieck basis via WC graph enumeration.

<a id="schubmult.rings.free_algebra.grove_basis.GroveBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from GroveBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.z_basis"></a>

# schubmult.rings.free\_algebra.z\_basis

<a id="schubmult.rings.free_algebra.z_basis.ZBasis"></a>

## ZBasis Objects

```python
class ZBasis(FreeAlgebraBasis)
```

Z basis of the free algebra.

Keys are tuples of positive integers (no zeros). The Z basis is
related to the Schubert basis by incrementing/decrementing code
entries by 1 and dropping zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a Z basis key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.pare_schubert"></a>

#### pare\_schubert

```python
@staticmethod
def pare_schubert(perm)
```

Extract the nonzero trimcode of *perm*, or None if it contains interior zeros.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Multiply two Z basis keys via shifted Schubert multiplication.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Z``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.z_basis.ZBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from ZBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.nelementary_basis"></a>

# schubmult.rings.free\_algebra.nelementary\_basis

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis"></a>

## NElementaryBasis Objects

```python
class NElementaryBasis(FreeAlgebraBasis)
```

Non-commutative elementary basis (L basis) of the free algebra.

Keys are tuples of positive integers. The transition to the word
basis uses composition refinements from SageMath.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.product"></a>

#### product

```python
@classmethod
def product(cls, key1, key2, coeff=S.One)
```

Concatenate two keys.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``L``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
def transition_word(cls, tup)
```

Transition an NElementary key to the word basis via composition refinements (requires SageMath).

<a id="schubmult.rings.free_algebra.nelementary_basis.NElementaryBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from NElementaryBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.jt_basis"></a>

# schubmult.rings.free\_algebra.jt\_basis

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis"></a>

## JTBasis Objects

```python
class JTBasis(FreeAlgebraBasis)
```

JT basis of the free algebra (J basis with a parameter *t*).

Keys are ``(tuple, int)`` pairs where the tuple is a nonzero code
and the integer tracks a power of the parameter *t*.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a JT key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.pare_schubert"></a>

#### pare\_schubert

```python
@staticmethod
def pare_schubert(perm)
```

Extract a nonzero trimcode from *perm*, or None if it contains zeros.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.normalize_dct"></a>

#### normalize\_dct

```python
@staticmethod
def normalize_dct(dct)
```

Normalize a word dict by collecting zeros into leading positions.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a *t*-weighted ``JT``-labelled display object.

<a id="schubmult.rings.free_algebra.jt_basis.JTBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from JTBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.j_basis"></a>

# schubmult.rings.free\_algebra.j\_basis

<a id="schubmult.rings.free_algebra.j_basis.JBasis"></a>

## JBasis Objects

```python
class JBasis(FreeAlgebraBasis)
```

J basis of the free algebra.

Keys are tuples of positive integers (no zeros allowed in transitions).
The J basis indexes elements whose Schubert expansion has no zero
entries in the Lehmer code.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.from_perm"></a>

#### from\_perm

```python
@staticmethod
def from_perm(perm, n)
```

Extract a J basis key from *perm* if the first *n* code entries are nonzero.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.coproduct"></a>

#### coproduct

```python
@classmethod
def coproduct(cls, key)
```

Coproduct for JBasis equals the bar-coproduct.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``J``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.j_basis.JBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from JBasis to *other_basis*.

<a id="schubmult.rings.free_algebra._core"></a>

# schubmult.rings.free\_algebra.\_core

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement"></a>

## FreeAlgebraElement Objects

```python
class FreeAlgebraElement(BaseRingElement)
```

Element of a free algebra, stored as a dict mapping basis keys to coefficients.

Keys are tuples of nonnegative integers (words in the word basis) or
basis-specific keys depending on the parent ring's basis. Supports
arithmetic operations, basis changes, and word-level operations like
injection, prefix, suffix, and interval extraction.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.interleave"></a>

#### interleave

```python
def interleave(other, zero_pad=True)
```

Interleave two elements letter-by-letter in the word basis.

Converts both elements to WordBasis, then interleaves each pair of
words by alternating entries (a1, b1, a2, b2, ...). Shorter words
are zero-padded when ``zero_pad`` is True.

**Arguments**:

- `other` - Another FreeAlgebraElement to interleave with.
- `zero_pad` - If True, pad shorter words with zeros.
  

**Returns**:

  The interleaved element in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.inject"></a>

#### inject

```python
def inject(i, other)
```

Insert another element's words at position *i* in this element's words.

Delegates to the current basis's ``inject`` classmethod.

**Arguments**:

- `i` - Nonnegative integer insertion index.
- `other` - Another FreeAlgebraElement to inject.
  

**Returns**:

  A new element with the injected words.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.prefix"></a>

#### prefix

```python
def prefix(length)
```

Extract the first *length* letters of each word.

Delegates to the current basis's ``prefix`` classmethod.

**Arguments**:

- `length` - Nonnegative integer prefix length.
  

**Returns**:

  A new element containing the prefixes.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.suffix"></a>

#### suffix

```python
def suffix(length)
```

Extract the last *length* letters of each word.

Delegates to the current basis's ``suffix`` classmethod.

**Arguments**:

- `length` - Nonnegative integer suffix length.
  

**Returns**:

  A new element containing the suffixes.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.interval"></a>

#### interval

```python
def interval(start, stop)
```

Extract a subword from position *start* to *stop* in each word.

Delegates to the current basis's ``interval`` classmethod.

**Arguments**:

- `start` - Nonnegative start index (inclusive).
- `stop` - Stop index (exclusive).
  

**Returns**:

  A new element containing the subwords.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.poly_inner_product"></a>

#### poly\_inner\_product

```python
def poly_inner_product(poly, genset, n)
```

Compute the inner product of this element with a polynomial.

Converts to WordBasis and pairs coefficient-by-coefficient with
the monomial expansion of *poly* in *genset*.

**Arguments**:

- `poly` - A polynomial expression.
- `genset` - The generating set of variables for the polynomial.
- `n` - Number of variables to use (or None for automatic).
  

**Returns**:

  The integer inner product value.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.kill_zero"></a>

#### kill\_zero

```python
def kill_zero(fat=False, val=S.Zero)
```

Remove zeros from each word key.

In the word basis, strips all zero entries from each key. When
*fat* is True, multiplies the coefficient by ``val`` raised to
the number of removed zeros instead of simply dropping them.

**Arguments**:

- `fat` - If True, weight by ``val`` per removed zero.
- `val` - The value to raise per zero when *fat* is True.
  

**Returns**:

  A new element with zeros removed, in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

Return a dict mapping printing terms to sympified coefficients.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.expand"></a>

#### expand

```python
def expand(deep=True, *args, **kwargs)
```

Expand all coefficients symbolically.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.pairing"></a>

#### pairing

```python
def pairing(other)
```

Compute the pairing of this element with *other* via the monomial basis.

Converts *self* to WordBasis and *other* to MonomialBasis, then
sums products of matching coefficients.

**Arguments**:

- `other` - Another element to pair with.
  

**Returns**:

  The integer pairing value.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.hom_nsym"></a>

#### hom\_nsym

```python
def hom_nsym()
```

Apply the homomorphism to noncommutative symmetric functions.

Strips zero entries from each word key (dropping them from the word)
and returns the result in the original basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_double_expand"></a>

#### tup\_double\_expand

```python
@staticmethod
@cache
def tup_double_expand(tup)
```

Expand a word tuple into the double Schubert basis via Pieri products.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.tup_expand"></a>

#### tup\_expand

```python
@staticmethod
@cache
def tup_expand(tup)
```

Expand a word tuple into the single Schubert basis via divide-and-conquer Pieri products.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis)
```

Convert this element to another basis.

**Arguments**:

- `other_basis` - The target basis class (e.g. WordBasis, SchubertBasis).
  

**Returns**:

  A new FreeAlgebraElement in the target basis's ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_expand"></a>

#### schub\_expand

```python
def schub_expand()
```

Expand this element into a single Schubert polynomial ring element.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.schub_double_expand"></a>

#### schub\_double\_expand

```python
def schub_double_expand()
```

Expand this element into a double Schubert polynomial ring element.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.bcoproduct"></a>

#### bcoproduct

```python
def bcoproduct()
```

Compute the bar-coproduct of this element in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.factorize"></a>

#### factorize

```python
def factorize(j)
```

Split each word at position *j*, returning a tensor element.

In the word basis, each word ``w`` maps to ``(w[:j], w[j:])``.
The result is expressed in the tensor ring of the original basis.

**Arguments**:

- `j` - Position at which to split each word.
  

**Returns**:

  An element of the tensor product ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.remove_zeros"></a>

#### remove\_zeros

```python
def remove_zeros(inserter=S.One)
```

Remove zero entries from each key, weighting by *inserter* per zero removed.

**Arguments**:

- `inserter` - Scalar multiplied per removed zero (default 1).
  

**Returns**:

  A new element with zeros stripped from keys.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.split"></a>

#### split

```python
def split(p)
```

Split each word at position *p* into a tensor element.

Words shorter than *p* are placed entirely in the left factor.

**Arguments**:

- `p` - Position at which to split.
  

**Returns**:

  An element of the tensor product ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebraElement.to_schub"></a>

#### to\_schub

```python
def to_schub(sym=False)
```

Convert to a Schubert ring element via ``tup_to_schub``.

**Arguments**:

- `sym` - If True, use symmetric expansion.
  

**Returns**:

  An element of the single Schubert ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra"></a>

## FreeAlgebra Objects

```python
class FreeAlgebra(BaseRing)
```

Free algebra ring with a configurable basis.

The algebra operates on :class:`FreeAlgebraElement` instances whose keys
are determined by the chosen basis (default :class:`WordBasis`). Supports
multiplication, tensor products, coproducts, and basis changes.

**Arguments**:

- `basis` - The basis class to use (default ``WordBasis``).
- `domain` - Coefficient domain (default ``EXRAW``).

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.mul_expr"></a>

#### mul\_expr

```python
def mul_expr(elem, x)
```

Multiply every coefficient of *elem* by the scalar *x*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.tensor_schub_expand"></a>

#### tensor\_schub\_expand

```python
def tensor_schub_expand(tensor)
```

Expand a tensor element into the Schubert polynomial tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.__init__"></a>

#### \_\_init\_\_

```python
def __init__(basis=WordBasis, domain=None)
```

Initialize a FreeAlgebra with the given basis and coefficient domain.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.right_pad"></a>

#### right\_pad

```python
@staticmethod
def right_pad(tup, n)
```

Right-pad *tup* with zeros to length *n*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(key)
```

Compute the coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.bcoproduct_on_basis"></a>

#### bcoproduct\_on\_basis

```python
@cache
def bcoproduct_on_basis(key)
```

Compute the bar-coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via the basis product rule.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.from_rc_graph"></a>

#### from\_rc\_graph

```python
def from_rc_graph(rc_graph)
```

Create an element from an RC graph.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.matmul"></a>

#### matmul

```python
def matmul(elem, other)
```

Internal product (``@`` operator) or scalar multiplication.

If *other* is a scalar, multiplies all coefficients. If *other* is
a FreeAlgebraElement, computes the internal product via the basis.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.new"></a>

#### new

```python
def new(*x)
```

Create a new element from the given key or arguments.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Return the display symbol for basis key *k*.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.from_dict"></a>

#### from\_dict

```python
def from_dict(element)
```

Construct an element from a dict of ``{key: coefficient}`` pairs.

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.skew_element"></a>

#### skew\_element

```python
def skew_element(w, u, n)
```

Skew schubert by elem sym

<a id="schubmult.rings.free_algebra._core.FreeAlgebra.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce a raw value into the coefficient domain.

<a id="schubmult.rings.free_algebra._core.make_AGx"></a>

#### make\_AGx

```python
def make_AGx(beta=None)
```

Create a Grothendieck-basis free algebra, optionally with custom beta.

<a id="schubmult.rings.free_algebra.lascoux_basis"></a>

# schubmult.rings.free\_algebra.lascoux\_basis

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis"></a>

## LascouxBasis Objects

```python
class LascouxBasis(FreeAlgebraBasis)
```

Lascoux polynomial (Demazure character) basis of the free algebra.

Lascouxs are weak composition tuples. Transitions to the Schubert basis
use RC graph enumeration filtered by extremal weight.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the LascouxPolyBasis as the dual of LascouxBasis.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Lascoux``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a Lascoux composition to the Grothendieck basis via WC graphs.

<a id="schubmult.rings.free_algebra.lascoux_basis.LascouxBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from LascouxBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.glide_basis"></a>

# schubmult.rings.free\_algebra.glide\_basis

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis"></a>

## GlideBasis Objects

```python
class GlideBasis(FreeAlgebraBasis)
```

Glide basis of the free algebra.

Keys are weak composition tuples. Transitions are computed via
polynomial algebra slide polynomial expansions.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``FS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
@classmethod
def transition_grothendieck(cls, key)
```

Transition a glide key to the Grothendieck basis.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the GlidePolyBasis as the dual.

<a id="schubmult.rings.free_algebra.glide_basis.GlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from GlideBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.key_basis"></a>

# schubmult.rings.free\_algebra.key\_basis

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis"></a>

## KeyBasis Objects

```python
class KeyBasis(FreeAlgebraBasis)
```

Key polynomial (Demazure character) basis of the free algebra.

Keys are weak composition tuples. Transitions to the Schubert basis
use RC graph enumeration filtered by extremal weight.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the KeyPolyBasis as the dual of KeyBasis.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return a ``Key``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a key composition to the Schubert basis via RC graphs.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.product"></a>

#### product

```python
@classmethod
@cache
def product(cls, key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to WordBasis and back.

<a id="schubmult.rings.free_algebra.key_basis.KeyBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from KeyBasis to *other_basis*.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis"></a>

# schubmult.rings.free\_algebra.fundamental\_slide\_basis

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis"></a>

## FundamentalSlideBasis Objects

```python
class FundamentalSlideBasis(FreeAlgebraBasis)
```

Fundamental slide basis of the free algebra.

Keys are weak composition tuples. Transitions are computed via
polynomial algebra slide polynomial expansions.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.is_key"></a>

#### is\_key

```python
@classmethod
def is_key(cls, x)
```

Return True if *x* is a tuple or list.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.as_key"></a>

#### as\_key

```python
@classmethod
def as_key(cls, x)
```

Normalize *x* to a tuple key.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.printing_term"></a>

#### printing\_term

```python
@classmethod
def printing_term(cls, k)
```

Return an ``FS``-labelled display object for key *k*.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition_schubert"></a>

#### transition\_schubert

```python
@classmethod
def transition_schubert(cls, key)
```

Transition a fundamental slide key to the Schubert basis.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition_word"></a>

#### transition\_word

```python
@classmethod
@cache
def transition_word(cls, key)
```

Transition a fundamental slide key to the word basis.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the FundamentalSlidePolyBasis as the dual.

<a id="schubmult.rings.free_algebra.fundamental_slide_basis.FundamentalSlideBasis.transition"></a>

#### transition

```python
@classmethod
def transition(cls, other_basis)
```

Return a transition function from FundamentalSlideBasis to *other_basis*.

<a id="schubmult.rings.schubert.double_schubert_ring"></a>

# schubmult.rings.schubert.double\_schubert\_ring

<a id="schubmult.rings.schubert.double_schubert_ring.DoubleSchubertElement"></a>

## DoubleSchubertElement Objects

```python
class DoubleSchubertElement(BaseSchubertElement)
```

Algebra with sympy coefficients
and a dict basis

<a id="schubmult.rings.schubert.double_grothendieck_ring"></a>

# schubmult.rings.schubert.double\_grothendieck\_ring

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckElement"></a>

## DoubleGrothendieckElement Objects

```python
class DoubleGrothendieckElement(BaseSchubertElement)
```

Element of a DoubleGrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing"></a>

## DoubleGrothendieckRing Objects

```python
class DoubleGrothendieckRing(BaseSchubertRing)
```

Ring of double (K-theoretic) Grothendieck polynomials G_w(x, y).

Elements are stored in the G-basis as ``{Permutation: coeff}``. There is
no direct structure-constant formula for the product implemented here:
instead both factors are expanded into the underlying ``DoubleSchubertRing``
(via ``grothendieck_poly_with_ring``), multiplied there, and the product is
converted back to the G-basis with ``to_groth_with_ring``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_root"></a>

#### exp\_root

```python
@cache
def exp_root(a, b)
```

Polynomial form of the root ``eps_a - eps_b``, i.e. ``(x^root - 1)/beta``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_weight"></a>

#### exp\_weight

```python
@cache
def exp_weight(index)
```

Polynomial form of the fundamental weight ``eps_index``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.exp_character"></a>

#### exp\_character

```python
def exp_character(weight)
```

Polynomial form of the multiplicative character ``x^weight``.

Determined by ``exp_root``: ``x^{eps_i} = 1 + beta * y_i``, so that
``(x^{eps_a - eps_b} - 1)/beta`` reproduces ``exp_root(a, b)``.

<a id="schubmult.rings.schubert.double_grothendieck_ring.DoubleGrothendieckRing.chevalley"></a>

#### chevalley

```python
def chevalley(weight, perm, n=None)
```

Lenart--Postnikov K_T-Chevalley formula for ``e^weight * G_perm``.

``weight`` is an integer vector in the ``epsilon`` basis.

<a id="schubmult.rings.schubert.chevalley"></a>

# schubmult.rings.schubert.chevalley

Lenart--Postnikov :math:`K_T`-Chevalley formula in type :math:`A`.

Implements Theorem 6.1 / Proposition 14.5 of Lenart--Postnikov, *Affine Weyl
groups in K-theory and representation theory* (arXiv:math/0309207)::

    e^lambda . [O_u] = sum_{w, mu} c^{lambda, mu}_{u, w} x^mu [O_w]

    c^{lambda, mu}_{u, w} = sum_J (-1)^{n(J)}

summed over subsets ``J = {j_1 < ... < j_s}`` of a fixed lambda-chain such that

(a) ``u > u r_{j_1} > ... > u r_{j_1} ... r_{j_s} = w`` is a saturated
    decreasing chain in Bruhat order, and
(b) ``-mu = u r_{j_1} ... r_{j_s} (-lambda)``,

with ``n(J)`` the number of negative roots among ``beta_{j_1}, ..., beta_{j_s}``.

Weights are integer vectors in the ``epsilon`` basis; a root ``(a, b)`` denotes
``eps_a - eps_b`` and is positive iff ``a < b``.

<a id="schubmult.rings.schubert.chevalley.fundamental_weight"></a>

#### fundamental\_weight

```python
def fundamental_weight(k, n)
```

``omega_k = eps_1 + ... + eps_k`` as a length-``n`` vector.

<a id="schubmult.rings.schubert.chevalley.lambda_chain"></a>

#### lambda\_chain

```python
@cache
def lambda_chain(weight, n)
```

Reduced ``lambda``-chain for ``weight`` in ``A_{n-1}``, via Prop. 6.7.

Returns a tuple of ``(beta, root, k)``: ``root`` is the positive root
``alpha`` of the affine reflection ``r_j = s_{alpha, k}``, and ``beta`` is
the signed root ``b(r_j)`` whose sign determines ``(-1)^{n(J)}``.

<a id="schubmult.rings.schubert.chevalley.kt_chevalley_coefficients"></a>

#### kt\_chevalley\_coefficients

```python
def kt_chevalley_coefficients(u, weight, n=None)
```

Chevalley coefficients for ``e^weight . [O_u]``.

Returns ``{(w, mu): coefficient}`` with ``mu`` an integer weight vector.

<a id="schubmult.rings.schubert.grothendieck_ring"></a>

# schubmult.rings.schubert.grothendieck\_ring

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckElement"></a>

## GrothendieckElement Objects

```python
class GrothendieckElement(BaseSchubertElement)
```

Element of a GrothendieckRing, stored as {Permutation: coeff}.

<a id="schubmult.rings.schubert.grothendieck_ring.GrothendieckRing"></a>

## GrothendieckRing Objects

```python
class GrothendieckRing(BaseSchubertRing)
```

Ring of Grothendieck polynomials.

A deformation of the Schubert polynomial ring with parameter beta.
Basis elements G_w satisfy G_u * G_v = sum_w c^w_{u,v}(beta) G_w
where c^w_{u,v}(beta) are polynomials in beta with integer coefficients.
When beta=0, recovers ordinary Schubert polynomials.

Parameters
----------
genset : GeneratingSet
    The generating set (variable alphabet).
beta : sympy/symengine symbol, optional
    The deformation parameter. Defaults to Symbol("β").

<a id="schubmult.rings.direct_product_ring"></a>

# schubmult.rings.direct\_product\_ring

<a id="schubmult.rings.direct_product_ring.DirectProductRing"></a>

## DirectProductRing Objects

```python
class DirectProductRing(BaseRing)
```

Direct product of an arbitrary (fixed) number of rings.

An element is stored as a dict mapping keys ``(i, k)`` to coefficients,
where ``i`` is the component index and ``k`` is a basis key from the
*i*-th ring.  Addition and multiplication are componentwise: terms from
different components never interact, and cross-component products are
zero.

Parameters
----------
``*rings`` : BaseRing
    One or more constituent rings.

Examples
--------
>>> D = DirectProductRing(R1, R2, R3)
>>> a = D.from_component(0, some_R1_element)
>>> b = D.from_component(2, some_R3_element)
>>> a + b          # lives in components 0 and 2
>>> a * b          # zero (different components)
>>> D[0]           # R1
>>> D.project(a, 0)  # back to R1

<a id="schubmult.rings.direct_product_ring.DirectProductRing.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(i)
```

Return the *i*-th constituent ring.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.component_one"></a>

#### component\_one

```python
def component_one(i)
```

Return the identity element supported only on component *i*.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.from_component"></a>

#### from\_component

```python
def from_component(i, elem)
```

Lift an element of ``self[i]`` into the direct product.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.project"></a>

#### project

```python
def project(elem, i)
```

Project onto component *i*, returning an element of ``self[i]``.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Componentwise multiplication.

Only terms in the same component interact; cross-component products
are zero.

<a id="schubmult.rings.direct_product_ring.DirectProductRing.__call__"></a>

#### \_\_call\_\_

```python
def __call__(*elems)
```

Construct an element from one element per component.

``D(e0, e1, ..., en)`` lifts each ``ei`` (an element of ``D[i]``)
into the direct product and sums them.

<a id="schubmult.rings.direct_product_ring.DirectProductRingElement"></a>

## DirectProductRingElement Objects

```python
class DirectProductRingElement(BaseRingElement)
```

<a id="schubmult.rings.direct_product_ring.DirectProductRingElement.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key)
```

Index by component integer or by ``(i, k)`` basis key.

* ``elem[i]`` — project onto component *i* (returns a ``self.ring[i]`` element).
* ``elem[(i, k)]`` — coefficient lookup (standard dict behaviour).

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.schubert\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis"></a>

## SchubertPolyBasis Objects

```python
class SchubertPolyBasis(PolynomialBasis)
```

Schubert polynomial basis.

Keys are ``(Permutation, length)`` pairs. Schubert polynomials form
the canonical basis for the polynomial algebra in Schubert calculus,
dual to the :class:`SchubertBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the coproduct of a Schubert key by splitting variable sets.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two Schubert keys using the Schubert ring multiplication.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_grothendieck"></a>

#### transition\_grothendieck

```python
def transition_grothendieck(dct)
```

Transition a Schubert dict to the Grothendieck polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_sepdesc"></a>

#### transition\_sepdesc

```python
def transition_sepdesc(dct, other_basis)
```

Transition from Schubert basis to separated descents basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_elementary"></a>

#### transition\_elementary

```python
def transition_elementary(dct, other_basis)
```

Transition from Schubert basis to elementary symmetric basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key_fundamental_slide"></a>

#### transition\_key\_fundamental\_slide

```python
def transition_key_fundamental_slide(perm, n)
```

Decompose a Schubert polynomial into fundamental slide polynomials via quasi-Yamanouchi RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
def transition_fundamental_slide(dct)
```

Transition a Schubert dict to the fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key_key"></a>

#### transition\_key\_key

```python
def transition_key_key(key)
```

Decompose a Schubert polynomial into key polynomials via highest-weight RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_key"></a>

#### transition\_key

```python
def transition_key(dct)
```

Transition a Schubert dict to the key polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a Schubert key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`SchubertBasis`).

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_forest_key"></a>

#### transition\_forest\_key

```python
def transition_forest_key(key)
```

Decompose a Schubert polynomial into forest polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition_forest"></a>

#### transition\_forest

```python
def transition_forest(dct)
```

Transition a Schubert dict to the forest polynomial basis.

<a id="schubmult.rings.polynomial_algebra.schubert_poly_basis.SchubertPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Schubert basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.composition_schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.composition\_schubert\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.composition_schubert_poly_basis.CompositionSchubertPolyBasis"></a>

## CompositionSchubertPolyBasis Objects

```python
class CompositionSchubertPolyBasis(SchubertPolyBasis)
```

Wrapper basis for Schubert polynomials indexed by weak compositions.

Keys are weak compositions interpreted as Lehmer codes. The underlying
computations are delegated to :class:`SchubertPolyBasis`, while display
uses the same Schubert printing as the corresponding permutation key.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.elem\_sym\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis"></a>

## ElemSymPolyBasis Objects

```python
class ElemSymPolyBasis(PolynomialBasis)
```

Elementary symmetric polynomial basis.

Keys are tuples encoding products of elementary symmetric polynomials
e_k(x_1, ..., x_n). Each key specifies degrees and variable counts
for the elementary symmetric factors.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition from elementary symmetric basis to Schubert basis.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from elementary symmetric basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.elem_sym_poly_basis.ElemSymPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from this basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.anti\_schubert\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis"></a>

## AntiSchubertPolyBasis Objects

```python
class AntiSchubertPolyBasis(PolynomialBasis)
```

Anti-Schubert polynomial basis.

Keys are ``(Permutation, length)`` pairs. This basis reverses the
monomial ordering relative to the standard Schubert basis, with
the coproduct correspondingly reversed.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the reversed coproduct of an anti-Schubert key.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two anti-Schubert keys using the underlying Schubert ring.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_key_key"></a>

#### transition\_key\_key

```python
def transition_key_key(key)
```

Decompose an anti-Schubert polynomial into key polynomials with reversed weights.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_key"></a>

#### transition\_key

```python
def transition_key(dct)
```

Transition an anti-Schubert dict to the key polynomial basis.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand an anti-Schubert key into reversed monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition_forest"></a>

#### transition\_forest

```python
def transition_forest(dct)
```

Transition an anti-Schubert dict to the forest polynomial basis.

<a id="schubmult.rings.polynomial_algebra.anti_schubert_poly_basis.AntiSchubertPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from anti-Schubert basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.lascoux\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis"></a>

## LascouxPolyBasis Objects

```python
class LascouxPolyBasis(PolynomialBasis)
```

Lascoux polynomial basis.

Keys are weak compositions. Lascoux polynomials provide a
basis that refines Grothendieck polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a glide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`LascouxBasis`).

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a Lascoux basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from Lascoux basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_glide_key"></a>

#### transition\_glide\_key

```python
def transition_glide_key(key)
```

Transition a Lascoux key to the glide basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition_glide"></a>

#### transition\_glide

```python
def transition_glide(dct)
```

Transition from Lascoux basis to glide basis.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Lascoux basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.lascoux_poly_basis.LascouxPolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two Lascoux keys using the Lascoux product rule.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.grothendieck\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis"></a>

## GrothendieckPolyBasis Objects

```python
class GrothendieckPolyBasis(PolynomialBasis)
```

Grothendieck polynomial basis.

Keys are ``(Permutation, length)`` pairs. Grothendieck polynomials form
the canonical basis for the polynomial algebra in Grothendieck calculus,
dual to the :class:`GrothendieckBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two Grothendieck keys using the Grothendieck ring multiplication.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition from Grothendieck basis to separated descents basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_glide_key"></a>

#### transition\_glide\_key

```python
def transition_glide_key(key)
```

Decompose a Grothendieck polynomial into glide polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_glide"></a>

#### transition\_glide

```python
def transition_glide(dct)
```

Transition a Grothendieck dict to the glide polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a Grothendieck key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GrothendieckBasis`).

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_grove_key"></a>

#### transition\_grove\_key

```python
def transition_grove_key(key)
```

Decompose a Grothendieck polynomial into grove polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_grove"></a>

#### transition\_grove

```python
def transition_grove(dct)
```

Transition a Grothendieck dict to the grove polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_lascoux_key"></a>

#### transition\_lascoux\_key

```python
def transition_lascoux_key(key)
```

Decompose a Grothendieck polynomial into Lascoux polynomials via omega insertion on RC-graphs.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition_lascoux"></a>

#### transition\_lascoux

```python
def transition_lascoux(dct)
```

Transition a Grothendieck dict to the Lascoux polynomial basis.

<a id="schubmult.rings.polynomial_algebra.grothendieck_poly_basis.GrothendieckPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from Grothendieck basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis"></a>

# schubmult.rings.polynomial\_algebra.base\_polynomial\_basis

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis"></a>

## PolynomialBasis Objects

```python
class PolynomialBasis(ABC)
```

Abstract base class for polynomial algebra bases.

Subclasses define how keys are represented, how to transition between
bases, and how to expand elements into explicit polynomials. Default
implementations delegate through the :class:`MonomialBasis`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.is_key"></a>

#### is\_key

```python
@abstractmethod
def is_key(x)
```

Return True if *x* is a valid key for this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.as_key"></a>

#### as\_key

```python
@abstractmethod
def as_key(x)
```

Normalize *x* into a canonical key for this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.attach_key"></a>

#### attach\_key

```python
def attach_key(dct)
```

Normalize all keys in *dct* via :meth:`as_key`.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.transition"></a>

#### transition

```python
@abstractmethod
def transition(other_basis)
```

Return a function mapping dicts of this basis to dicts in *other_basis*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.from_expr"></a>

#### from\_expr

```python
def from_expr(expr, length=None)
```

Parse a symbolic expression into this basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.printing_term"></a>

#### printing\_term

```python
@abstractmethod
def printing_term(k)
```

Return the display symbol for key *k*.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.compose_transition"></a>

#### compose\_transition

```python
@staticmethod
def compose_transition(tkeyfunc, output)
```

Apply a transition function to a dict of basis elements.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.change_tensor_basis"></a>

#### change\_tensor\_basis

```python
@classmethod
def change_tensor_basis(cls, tensor_elem, basis1, basis2)
```

Change the bases of both factors of a tensor element.

**Arguments**:

- `tensor_elem` - An element of a tensor product ring.
- `basis1` - Target basis for the left factor.
- `basis2` - Target basis for the right factor.
  

**Returns**:

  The tensor element re-expressed in the new bases.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a basis dict into an explicit polynomial expression.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the coproduct of *key* by delegating through the monomial basis.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two keys by transitioning to the monomial basis and back.

<a id="schubmult.rings.polynomial_algebra.base_polynomial_basis.PolynomialBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.fundamental\_slide\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.get_descent_composition"></a>

#### get\_descent\_composition

```python
def get_descent_composition(word)
```

Compute the descent composition of a word.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.slide_product"></a>

#### slide\_product

```python
def slide_product(a, b)
```

Compute the structure constants for multiplying two fundamental slide polynomials.

Given weak compositions *a* and *b*, returns a dict mapping result
compositions to their coefficients in the fundamental slide expansion
of the product.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis"></a>

## FundamentalSlidePolyBasis Objects

```python
class FundamentalSlidePolyBasis(PolynomialBasis)
```

Fundamental slide polynomial basis.

Keys are weak compositions. Fundamental slide polynomials provide a
basis that refines Schubert polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a slide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`FundamentalSlideBasis`).

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a slide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from fundamental slide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from fundamental slide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis.FundamentalSlidePolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two fundamental slide keys using the slide product rule.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.double\_forest\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis"></a>

## DoubleForestPolyBasis Objects

```python
class DoubleForestPolyBasis(PolynomialBasis)
```

Abstract double forest polynomial basis.

Keys are weak compositions indexing double forest basis elements DF[key],
with polynomial coefficients in a second generating set (equivariant vars).

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis.basis_polynomial"></a>

#### basis\_polynomial

```python
def basis_polynomial(key)
```

Expand one double-forest basis key as a polynomial in x,t.

<a id="schubmult.rings.polynomial_algebra.double_forest_poly_basis.DoubleForestPolyBasis.basis_forest_expansion"></a>

#### basis\_forest\_expansion

```python
def basis_forest_expansion(key, length)
```

Expand one double-forest basis key into ForestPolyBasis in x.

<a id="schubmult.rings.polynomial_algebra._core"></a>

# schubmult.rings.polynomial\_algebra.\_core

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement"></a>

## PolynomialAlgebraElement Objects

```python
class PolynomialAlgebraElement(BaseRingElement)
```

Element of a polynomial algebra, stored as a dict mapping basis keys to coefficients.

Keys are exponent tuples (in the monomial basis) or basis-specific keys
depending on the parent ring's basis. Supports arithmetic, basis changes,
and duality pairing with free algebra elements.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.as_coefficients_dict"></a>

#### as\_coefficients\_dict

```python
def as_coefficients_dict()
```

Return a dict mapping printing terms to sympified coefficients.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.change_basis"></a>

#### change\_basis

```python
def change_basis(other_basis: type)
```

Convert this element to another polynomial basis.

**Arguments**:

- `other_basis` - A basis class, basis instance, or callable returning a basis.
  

**Returns**:

  A new PolynomialAlgebraElement in the target basis's ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.expand"></a>

#### expand

```python
def expand()
```

Expand this element into an explicit polynomial expression.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebraElement.apply_dual_element"></a>

#### apply\_dual\_element

```python
def apply_dual_element(dual_elem)
```

Pair this polynomial element with a dual free algebra element.

Converts *self* to the monomial basis and *dual_elem* to the word
basis, then sums products of matching coefficients.

**Arguments**:

- `dual_elem` - A FreeAlgebraElement to pair with.
  

**Returns**:

  The scalar pairing value.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra"></a>

## PolynomialAlgebra Objects

```python
class PolynomialAlgebra(BaseRing)
```

Polynomial algebra ring with a configurable basis.

The algebra operates on :class:`PolynomialAlgebraElement` instances whose
keys are determined by the chosen basis. Supports multiplication, basis
changes, coproducts, and conversion from symbolic expressions.

**Arguments**:

- `basis` - A basis instance (e.g. ``MonomialBasis(x)``).
- `domain` - Coefficient domain (default ``EXRAW``).

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.__init__"></a>

#### \_\_init\_\_

```python
def __init__(basis, domain=None)
```

Initialize a PolynomialAlgebra with the given basis and coefficient domain.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
@cache
def coproduct_on_basis(key)
```

Compute the coproduct of a single basis key in the tensor ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.mul"></a>

#### mul

```python
def mul(elem, other)
```

Multiply two elements via the basis product rule.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.new"></a>

#### new

```python
def new(*x)
```

Create a new element from the given key or expression.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.from_expr"></a>

#### from\_expr

```python
def from_expr(x, length=None)
```

Create an element from a symbolic expression.

Parses *x* into monomials, then transitions to this ring's basis.

**Arguments**:

- `x` - A symbolic polynomial expression.
- `length` - Optional fixed number of variables.
  

**Returns**:

  A PolynomialAlgebraElement in this ring.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.printing_term"></a>

#### printing\_term

```python
def printing_term(k)
```

Return the display symbol for basis key *k*.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.from_dict"></a>

#### from\_dict

```python
def from_dict(element)
```

Construct an element from a dict of ``{key: coefficient}`` pairs.

<a id="schubmult.rings.polynomial_algebra._core.PolynomialAlgebra.domain_new"></a>

#### domain\_new

```python
def domain_new(element, orig_domain=None)
```

Coerce a raw value into the coefficient domain.

<a id="schubmult.rings.polynomial_algebra.monomial_basis"></a>

# schubmult.rings.polynomial\_algebra.monomial\_basis

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis"></a>

## MonomialBasis Objects

```python
class MonomialBasis(PolynomialBasis)
```

Standard monomial basis for the polynomial algebra.

Keys are tuples of nonnegative integers representing exponent vectors.
This is the fundamental basis through which other bases transition
by default, and is dual to the :class:`WordBasis` of the free algebra.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.coproduct"></a>

#### coproduct

```python
def coproduct(key)
```

Compute the deconcatenation coproduct on a monomial key.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two monomial keys by component-wise addition of exponents.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.expand_monom"></a>

#### expand\_monom

```python
def expand_monom(monom)
```

Convert an exponent tuple to a monomial expression in the generating set.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a dict of monomial keys into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_slide"></a>

#### transition\_slide

```python
def transition_slide(dct, other_basis)
```

Transition a monomial dict to a slide-type basis via triangular inversion.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_slide_monom"></a>

#### transition\_slide\_monom

```python
def transition_slide_monom(other_basis, monom, coeff=S.One)
```

Express a single monomial in a slide-type basis by dominance-order inversion.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_anti_schubert"></a>

#### transition\_anti\_schubert

```python
def transition_anti_schubert(dct, other_basis)
```

Transition monomials to the anti-Schubert basis by reversing exponents.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct)
```

Transition monomials to the Schubert basis by grouping by length.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition_double_forest"></a>

#### transition\_double\_forest

```python
def transition_double_forest(dct, other_basis)
```

Transition monomials to DoubleForestPolyBasis via ForestPolyBasis.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`WordBasis`).

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from monomial basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.monomial_basis.MonomialBasis.from_expr"></a>

#### from\_expr

```python
def from_expr(expr, length=None)
```

Parse a symbolic expression into monomial-basis coefficient dict.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.key\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.traverse_demaz"></a>

#### traverse\_demaz

```python
def traverse_demaz(pl, w)
```

Traverse the Demazure graph starting from plactic element *pl* along the code word of *w*.

Yields all distinct elements reachable by successive lowering operators.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis"></a>

## KeyPolyBasis Objects

```python
class KeyPolyBasis(PolynomialBasis)
```

Key polynomial (Demazure character) basis.

Keys are weak compositions. Key polynomials are characters of Demazure
modules, computed by applying Demazure operators to highest-weight
plactic tableaux.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`KeyBasis`).

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a key polynomial key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a key basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from key basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.key_poly_basis.KeyPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from key basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.forest\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis"></a>

## ForestPolyBasis Objects

```python
class ForestPolyBasis(PolynomialBasis)
```

Forest polynomial basis.

Keys are weak compositions encoding indexed forests. Forest polynomials
are computed by summing over decreasing labelings of the corresponding
forest structure.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a forest key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a forest basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from forest basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition_fundamental_slide"></a>

#### transition\_fundamental\_slide

```python
def transition_fundamental_slide(dct)
```

Transition from forest basis to fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.to_fundamental_slide"></a>

#### to\_fundamental\_slide

```python
def to_fundamental_slide(key)
```

Express a single forest key in the fundamental slide basis.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`ForestBasis`).

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from forest basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.forest_poly_basis.ForestPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two forest keys by transitioning through the Schubert basis.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.monomial\_slide\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis"></a>

## MonomialSlidePolyBasis Objects

```python
class MonomialSlidePolyBasis(PolynomialBasis)
```

Monomial slide polynomial basis.

Keys are weak compositions. Monomial slide polynomials refine
key polynomials and are coarser than monomials, using a recursive
construction based on the first nonzero entry.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a monomial slide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from monomial slide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a monomial slide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.monomial_slide_poly_basis.MonomialSlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from monomial slide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.grove\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis"></a>

## GrovePolyBasis Objects

```python
class GrovePolyBasis(PolynomialBasis)
```

Grove polynomial basis.

Keys are weak compositions encoding indexed forests. Grove polynomials are
the ``beta``-deformed (set-valued) forest polynomials, computed by summing
over set-valued labelings of the corresponding forest structure.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a grove key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a grove basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from grove basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two grove keys using the glide product rule.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from grove basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.grove_poly_basis.GrovePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GroveBasis`).

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.glide\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.glide_monomials"></a>

#### glide\_monomials

```python
def glide_monomials(key)
```

Monomial expansion of the glide polynomial :math:`\mathcal{G}_{key}`.

Returns a dict mapping each exponent tuple ``v`` (a weak composition of the
same length as ``key``) to the integer coefficient of the monomial
:math:`x^v`, i.e. the number of glides of ``key`` with weight ``v``. The
corresponding power of ``beta`` for the weight ``v`` is
``sum(v) - sum(key)`` (the excess), which is constant across all glides of a
given weight, so it need not be stored explicitly.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.glide_product"></a>

#### glide\_product

```python
def glide_product(key1, key2)
```

Structure constants for a product of two glide polynomials.

Implements the Littlewood-Richardson rule of O. Pechenik and D. Searles,
"Decompositions of Grothendieck Polynomials" (arXiv:1611.02545), Theorem
4.9, which expands the product of the glide polynomials indexed by the weak
compositions ``key1`` and ``key2`` in the glide basis:

.. math::

    \mathcal{G}_a \, \mathcal{G}_b
        = \sum_c \beta^{|c| - |a| - |b|} \, g_{a,b}^{c} \, \mathcal{G}_c .

Rather than enumerating the genomic shuffle set directly, we compute the
(uniquely determined) coefficients by expanding the product in monomials and
straightening into the glide basis with the leading-term algorithm from the
proof that the glide polynomials form a basis (Theorem 2.6). Because the
excess of a glide of ``v`` equals ``sum(v) - sum(index)``, the power of
``beta`` is recovered from the total degree and only the positive integer
multiplicities :math:`g_{a,b}^{c}` are returned.

Both compositions are padded with trailing zeros to a common length ``n``;
every key ``c`` in the returned dict is a weak composition of length ``n``.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis"></a>

## GlidePolyBasis Objects

```python
class GlidePolyBasis(PolynomialBasis)
```

Glide polynomial basis.

Keys are weak compositions. Glide polynomials provide a
basis that refines Grothendieck polynomials and coarsens monomials, with
an efficient combinatorial product rule.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.to_monoms"></a>

#### to\_monoms

```python
def to_monoms(key)
```

Expand a glide key into a dict of monomial exponent tuples.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.dual_basis"></a>

#### dual\_basis

```python
@classmethod
def dual_basis(cls)
```

Return the dual free algebra basis class (:class:`GlideBasis`).

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.expand"></a>

#### expand

```python
def expand(dct)
```

Expand a glide basis dict into a symbolic polynomial expression.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.transition_monomial"></a>

#### transition\_monomial

```python
def transition_monomial(dct)
```

Transition from glide basis to monomial basis.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from glide basis to *other_basis*.

<a id="schubmult.rings.polynomial_algebra.glide_poly_basis.GlidePolyBasis.product"></a>

#### product

```python
@cache
def product(key1, key2, coeff=S.One)
```

Multiply two glide keys using the glide product rule.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis"></a>

# schubmult.rings.polynomial\_algebra.sepdesc\_poly\_basis

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis"></a>

## SepDescPolyBasis Objects

```python
class SepDescPolyBasis(PolynomialBasis)
```

Separated descents polynomial basis.

Keys are ``(Permutation, Permutation, k)`` triples. Elements are
products of pairs of Schubert polynomials, parameterized by a
separation level *k*.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.product"></a>

#### product

```python
def product(key1, key2, coeff=S.One)
```

Multiply two separated-descents keys by transitioning through Schubert.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.transition_schubert"></a>

#### transition\_schubert

```python
def transition_schubert(dct, other_basis)
```

Transition from separated descents to Schubert basis by multiplying the pair factors.

<a id="schubmult.rings.polynomial_algebra.sepdesc_poly_basis.SepDescPolyBasis.transition"></a>

#### transition

```python
def transition(other_basis)
```

Return a transition function from separated descents to *other_basis*.

<a id="schubmult.rings.printing"></a>

# schubmult.rings.printing

<a id="schubmult.rings.printing.SepDescSchubPoly"></a>

## SepDescSchubPoly Objects

```python
class SepDescSchubPoly(PrintingTerm)
```

Printing term for SeparatedDescentsRing: Xi_{perm}^{length}

<a id="schubmult.rings.printing.GrothendieckPoly"></a>

## GrothendieckPoly Objects

```python
class GrothendieckPoly(PrintingTerm)
```

Printing term for GrothendieckRing: G_w(x)

<a id="schubmult.rings.printing.DoubleGrothendieckPoly"></a>

## DoubleGrothendieckPoly Objects

```python
class DoubleGrothendieckPoly(PrintingTerm)
```

Printing term for DoubleGrothendieckRing: G_w(x, y)

<a id="schubmult.rings.base_ring"></a>

# schubmult.rings.base\_ring

<a id="schubmult.rings.base_ring.BaseRing"></a>

## BaseRing Objects

```python
class BaseRing(Ring, CompositeDomain)
```

<a id="schubmult.rings.base_ring.BaseRing.from_dict_unchecked"></a>

#### from\_dict\_unchecked

```python
def from_dict_unchecked(element)
```

from_dict for coefficients already known to lie in the domain (drops structural zeros only).

<a id="schubmult.rings.tensor_ring"></a>

# schubmult.rings.tensor\_ring

<a id="schubmult.rings.tensor_ring.TensorRing"></a>

## TensorRing Objects

```python
class TensorRing(BaseRing)
```

<a id="schubmult.rings.tensor_ring.TensorRing.coproduct_on_basis"></a>

#### coproduct\_on\_basis

```python
def coproduct_on_basis(k)
```

Compute coproduct of a basis element in the tensor ring.

For k = (k_1, k_2, ..., k_n) in ring_1 ⊗ ring_2 ⊗ ... ⊗ ring_n,
Δ(k_1 ⊗ k_2 ⊗ ... ⊗ k_n) = (⊗ Δ(k_i))

This properly interlaces the individual coproducts.

<a id="schubmult.rings.tensor_ring.TensorRing.mul"></a>

#### mul

```python
def mul(elem1, elem2)
```

Multiply two elements in the tensor ring.

(a1 ⊗ a2 ⊗ ... ⊗ an) * (b1 ⊗ b2 ⊗ ... ⊗ bn) = (a1*b1) ⊗ (a2*b2) ⊗ ... ⊗ (an*bn)

<a id="schubmult.rings.tensor_ring.TensorRingElement"></a>

## TensorRingElement Objects

```python
class TensorRingElement(BaseRingElement)
```

<a id="schubmult.rings.tensor_ring.TensorRingElement.coproduct"></a>

#### coproduct

```python
def coproduct()
```

Override coproduct to use the correct target ring.

<a id="schubmult.rings.combinatorial.plactic_algebra"></a>

# schubmult.rings.combinatorial.plactic\_algebra

<a id="schubmult.rings.combinatorial.plactic_algebra.PlacticAlgebraElement"></a>

## PlacticAlgebraElement Objects

```python
class PlacticAlgebraElement(BaseRingElement)
```

PlacticAlgebra elements are linear combinations of Plactic basis elements.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring"></a>

# schubmult.rings.combinatorial.schubert\_monomial\_ring

Schubert Monomial Ring module

Provides base classes for rings whose basis elements represent Schubert monomials
(e.g., RC-graphs, BPDs, pipe dreams) with common operations like expansion to
polynomials, divided differences, and crystal operations.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialPrintingTerm"></a>

## SchubertMonomialPrintingTerm Objects

```python
class SchubertMonomialPrintingTerm(TypedPrintingTerm)
```

Printing term for Schubert monomial basis elements.

Delegates printing to the underlying key object (typically an RCGraph, BPD, etc.)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement"></a>

## SchubertMonomialRingElement Objects

```python
class SchubertMonomialRingElement(BaseRingElement)
```

Base class for ring elements whose basis elements are Schubert monomials.

This provides a common interface for objects like:
- RCGraphRingElement (basis elements are RCGraphs)
- BPDRingElement (basis elements are BPDs)

Common operations include:
- Polynomial expansion via polyvalue()
- Divided difference operators
- Crystal structure operations (if the basis elements support them)

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x, y=None, *args, **kwargs)
```

Evaluate as a polynomial in variables x (and optionally y).

Linear extension: for each basis element, call its polyvalue() method
and sum the results weighted by coefficients.

**Arguments**:

- `x` - Variable or sequence of variables for polynomial evaluation
- `y` - Optional second set of variables for double Schubert polynomials
- ```**kwargs``` - Additional arguments passed to basis element polyvalue
  

**Returns**:

  Symbolic expression representing the polynomial

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRingElement.to_free_algebra_element"></a>

#### to\_free\_algebra\_element

```python
def to_free_algebra_element(basis=None, *, word=False)
```

Convert to FreeAlgebra element in Schubert basis.

<a id="schubmult.rings.combinatorial.schubert_monomial_ring.SchubertMonomialRing"></a>

## SchubertMonomialRing Objects

```python
class SchubertMonomialRing(BaseRing)
```

Base class for rings whose basis elements are Schubert monomials.

Inherits from BaseRing to provide standard ring operations (add, sub, mul, etc.)

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_rc\_factor\_algebra

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebraElement"></a>

## BoundedRCFactorAlgebraElement Objects

```python
class BoundedRCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebraElement.prune"></a>

#### prune

```python
def prune()
```

Merge terms whose evaluated RC graphs share forest/omega invariants.

Terms are grouped by
``(rc.forest_weight, rc.omega_invariant[1])`` where
``rc = self.ring.key_to_rc_graph(key)``.
One key per group is kept (first encountered), and coefficients are summed.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebra"></a>

## BoundedRCFactorAlgebra Objects

```python
class BoundedRCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.bounded_rc_factor_algebra.BoundedRCFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring"></a>

# schubmult.rings.combinatorial.dual\_rc\_graph\_ring

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement"></a>

## DualRCGraphRingElement Objects

```python
class DualRCGraphRingElement(SchubertMonomialRingElement)
```

DualRCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRing"></a>

## DualRCGraphRing Objects

```python
class DualRCGraphRing(SchubertMonomialRing)
```

<a id="schubmult.rings.combinatorial.dual_rc_graph_ring.DualRCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the DualRCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.wc_graph_ring"></a>

# schubmult.rings.combinatorial.wc\_graph\_ring

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRingElement"></a>

## WCGraphRingElement Objects

```python
class WCGraphRingElement(SchubertMonomialRingElement)
```

WCGraphRing elements are linear combinations of WCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing"></a>

## WCGraphRing Objects

```python
class WCGraphRing(SchubertMonomialRing)
```

<a id="schubmult.rings.combinatorial.wc_graph_ring.WCGraphRing.groth"></a>

#### groth

```python
def groth(perm, n=None)
```

Return the WCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.rc_schubert_ring"></a>

# schubmult.rings.combinatorial.rc\_schubert\_ring

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebraElement"></a>

## BoundedRCFactorAlgebraElement Objects

```python
class BoundedRCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebraElement.prune"></a>

#### prune

```python
def prune()
```

Merge terms whose evaluated RC graphs share forest/omega invariants.

Terms are grouped by
``(rc.forest_weight, rc.omega_invariant[1])`` where
``rc = self.ring.key_to_rc_graph(key)``.
One key per group is kept (first encountered), and coefficients are summed.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra"></a>

## BoundedRCFactorAlgebra Objects

```python
class BoundedRCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.rc_schubert_ring.BoundedRCFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.rc_graph_ring"></a>

# schubmult.rings.combinatorial.rc\_graph\_ring

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement"></a>

## RCGraphRingElement Objects

```python
class RCGraphRingElement(CrystalGraphRingElement, SchubertMonomialRingElement)
```

RCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.squash_product"></a>

#### squash\_product

```python
def squash_product(other_rc)
```

Linear extension of RCGraph.squash_product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.vertical_coproduct"></a>

#### vertical\_coproduct

```python
def vertical_coproduct()
```

Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
and induces the mul product.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an RCGraphRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an RCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an RCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRingElement.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

Linear extension of RCGraph.crystal_reflection:
For each basis RCGraph, apply its crystal_reflection(index) and collect results.

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing"></a>

## RCGraphRing Objects

```python
class RCGraphRing(SchubertMonomialRing, CrystalGraphRing)
```

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the RCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.rc_graph_ring.RCGraphRing.double_schub"></a>

#### double\_schub

```python
def double_schub(perm, coeff_genset, n=None)
```

Return the RCGraphRing element corresponding to the double Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used), with `coeff_genset`
as the set of variables for the double part.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring"></a>

# schubmult.rings.combinatorial.alt\_rc\_graph\_ring

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement"></a>

## AltRCGraphRingElement Objects

```python
class AltRCGraphRingElement(CrystalGraphRingElement,
                            SchubertMonomialRingElement)
```

AltRCGraphRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.divdiff_perm"></a>

#### divdiff\_perm

```python
def divdiff_perm(perm)
```

Apply divided difference operator for `perm` to self.
Linear extension of RCGraph.divdiff_perm.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.divdiff"></a>

#### divdiff

```python
def divdiff(*seq)
```

Sequential divided difference operators.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.vertical_coproduct"></a>

#### vertical\_coproduct

```python
def vertical_coproduct()
```

Coproduct of RC graphs, coincides with the coproduct on Schubert polynomials
and induces the mul product.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an AltRCGraphRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an AltRCGraphRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRingElement.crystal_reflection"></a>

#### crystal\_reflection

```python
def crystal_reflection(index)
```

Linear extension of RCGraph.crystal_reflection:
For each basis RCGraph, apply its crystal_reflection(index) and collect results.

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRing"></a>

## AltRCGraphRing Objects

```python
class AltRCGraphRing(SchubertMonomialRing, CrystalGraphRing)
```

<a id="schubmult.rings.combinatorial.alt_rc_graph_ring.AltRCGraphRing.schub"></a>

#### schub

```python
def schub(perm, n=None)
```

Return the AltRCGraphRing element corresponding to the Schubert polynomial
indexed by `perm` in `S_n` (if n is None, n = len(perm) is used).

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_rc\_forest\_factor\_algebra

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebraElement"></a>

## BoundedRCForestFactorAlgebraElement Objects

```python
class BoundedRCForestFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedRCForestFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra"></a>

## BoundedRCForestFactorAlgebra Objects

```python
class BoundedRCForestFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.bounded_rc_forest_factor_algebra.BoundedRCForestFactorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.crystal_graph_ring"></a>

# schubmult.rings.combinatorial.crystal\_graph\_ring

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRing"></a>

## CrystalGraphRing Objects

```python
class CrystalGraphRing(BaseRing)
```

Ring whose basis elements are CrystalGraph-like objects.

We deliberately do not special-case tensor objects here: CrystalGraphTensor
implements the same CrystalGraph API and will be handled by polymorphism.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement"></a>

## CrystalGraphRingElement Objects

```python
class CrystalGraphRingElement(BaseRingElement, CrystalGraph)
```

Element of the CrystalGraphRing.

Keys are arbitrary objects that implement the CrystalGraph API (including
CrystalGraphTensor). All crystal operators / statistics are lifted linearly
by delegating to the underlying key's methods.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index: int)
```

Linearized raising operator: delegate to each key's raising_operator
and collect results in the ring.

<a id="schubmult.rings.combinatorial.crystal_graph_ring.CrystalGraphRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index: int)
```

Linearized lowering operator: delegate to each key's lowering_operator.

<a id="schubmult.rings.combinatorial.eg_plactic_ring"></a>

# schubmult.rings.combinatorial.eg\_plactic\_ring

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement"></a>

## EGPlacticRingElement Objects

```python
class EGPlacticRingElement(CrystalGraphRingElement)
```

EGPlacticRing elements are linear combinations of RCGraph basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.__mod__"></a>

#### \_\_mod\_\_

```python
def __mod__(other)
```

Polynomial product: self % other.
Currently only defined when `other` is a dominant RC graph.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

Linear extension of RCGraph.raising_operator:
Apply raising_operator(index) to every basis RCGraph in self, collect results.
Returns an EGPlacticRingElement (possibly zero).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

Linear extension of RCGraph.lowering_operator.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.phi"></a>

#### phi

```python
def phi(index)
```

phi(element) := max_{basis rc in supp(element)} phi(rc)
If element is zero, returns 0.

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.epsilon"></a>

#### epsilon

```python
def epsilon(index)
```

epsilon(element) := max_{basis rc in supp(element)} epsilon(rc)

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Use maximum crystal length of basis graphs in support (0 for the zero element).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an EGPlacticRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight()
```

Iteratively raise the element until no further raising is possible.
Returns (highest_weight_element, raise_seq).

Behavior notes:
- This is the natural linear-extension of CrystalGraph.to_highest_weight.
- The returned `highest_weight_element` is an EGPlacticRingElement.
- raise_seq is the sequence of row indices applied (in order).

<a id="schubmult.rings.combinatorial.eg_plactic_ring.EGPlacticRingElement.reverse_raise_seq"></a>

#### reverse\_raise\_seq

```python
def reverse_raise_seq(raise_seq)
```

Apply lowering_operator in reverse order to `raise_seq`.
If the path dies (result is zero), return None (mirrors scalar behavior).

<a id="schubmult.rings.combinatorial.grass_tensor_algebra"></a>

# schubmult.rings.combinatorial.grass\_tensor\_algebra

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebraElement"></a>

## GrassTensorAlgebraElement Objects

```python
class GrassTensorAlgebraElement(CrystalGraphRingElement)
```

Element of GrassTensorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra"></a>

## GrassTensorAlgebra Objects

```python
class GrassTensorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian RC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
RC graph. Simplification rules

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra.dual_product_on_basis"></a>

#### dual\_product\_on\_basis

```python
def dual_product_on_basis(left_key, right_key)
```

Dual product to the coproduct_on_basis deconcatenation.

<a id="schubmult.rings.combinatorial.grass_tensor_algebra.GrassTensorAlgebra.key_to_rc_graph"></a>

#### key\_to\_rc\_graph

```python
def key_to_rc_graph(key: CrystalGraphTensor | tuple) -> RCGraph
```

Evaluate a tensor key to an RCGraph using left-to-right squash_product.

<a id="schubmult.rings.combinatorial.chute_move_ring"></a>

# schubmult.rings.combinatorial.chute\_move\_ring

<a id="schubmult.rings.combinatorial.chute_move_ring.ChuteMoveRingElement"></a>

## ChuteMoveRingElement Objects

```python
class ChuteMoveRingElement(SchubertMonomialRingElement)
```

ChuteMoveRing elements are linear combinations of ChuteMoveElement basis elements.

The product % is the polynomial product. Currently only defined when the right side
is a dominant RC graph.

The Leibniz rule should hold for % somehow. Claude's idea is to define the ambiguous term in the Leibniz formula instead of trying
to do this directly.

The product * is well defined for any pair of RC graphs and is the dual product.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra"></a>

# schubmult.rings.combinatorial.bounded\_wc\_factor\_algebra

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebraElement"></a>

## BoundedWCFactorAlgebraElement Objects

```python
class BoundedWCFactorAlgebraElement(CrystalGraphRingElement)
```

Element of BoundedWCFactorAlgebra: finite linear combinations of Grass tensors.

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebra"></a>

## BoundedWCFactorAlgebra Objects

```python
class BoundedWCFactorAlgebra(CrystalGraphRing)
```

Tensor-like algebra on tuples of full Grassmannian WC graphs.

Basis keys are tuples (g1, ..., gk) where each gi is a full Grassmannian
WC graph. Simplification rules

<a id="schubmult.rings.combinatorial.bounded_wc_factor_algebra.BoundedWCFactorAlgebra.key_to_wc_graph"></a>

#### key\_to\_wc\_graph

```python
def key_to_wc_graph(key) -> WCGraph
```

Evaluate a tensor key to an WCGraph using left-to-right squash_product.

<a id="schubmult.combinatorics.set_valued_tableau"></a>

# schubmult.combinatorics.set\_valued\_tableau

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau"></a>

## SetValuedTableau Objects

```python
class SetValuedTableau(GridPrint, CrystalGraph)
```

A semistandard set-valued tableau.

Each box of a Young diagram holds a non-empty *set* of positive integers.
The tableau is semistandard in the set-valued sense: reading each box as its
set of labels, ``max`` of a box is strictly less than ``min`` of the box
below it (columns strictly increase) and weakly less than ``min`` of the box
to its right (rows weakly increase).

Internally the tableau is stored as a dict ``{(row, col): tuple(labels)}``
where each ``labels`` tuple is sorted ascending. Empty boxes are simply
absent from the dict.

The tableau crystal operators are realized via the tensor model on
:class:`~schubmult.combinatorics.set_word.SetWord` with
:class:`~schubmult.combinatorics.set_word.SetLetter` factors, using
row-reading order over boxes (bottom-to-top, left-to-right).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.__init__"></a>

#### \_\_init\_\_

```python
def __init__(cells=None)
```

Create a set-valued tableau.

``cells`` may be:

- a dict ``{(row, col): iterable of labels}``, or
- a nested sequence of rows, where each entry is an iterable of labels
  (or a single integer for a singleton box); ``None`` marks an empty
  leading (skew) box.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.from_cells"></a>

#### from\_cells

```python
@classmethod
def from_cells(cls, cells)
```

Build a :class:`SetValuedTableau` from a dict ``{(row, col): labels}``.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.cells"></a>

#### cells

```python
@property
def cells()
```

Return the underlying ``{(row, col): tuple(labels)}`` dict (a copy).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.weight"></a>

#### weight

```python
@property
def weight()
```

Return the content weight: ``weight[v - 1]`` counts occurrences of ``v``.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.add_label"></a>

#### add\_label

```python
def add_label(row, col, label)
```

Return a new tableau with ``label`` added to the box at ``(row, col)``.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.is_semistandard"></a>

#### is\_semistandard

```python
def is_semistandard()
```

Check the semistandard set-valued conditions on rows and columns.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i)
```

Crystal lowering operator ``f_i`` via :class:`SetWord`.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

Crystal raising operator ``e_i`` via :class:`SetWord`.

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Content weight ``(`1`, `2`, ...)`` (alias of :attr:`weight`).

<a id="schubmult.combinatorics.set_valued_tableau.SetValuedTableau.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Upper bound on crystal operator indices (matches ``Plactic``).

<a id="schubmult.combinatorics.anti_rc_graph"></a>

# schubmult.combinatorics.anti\_rc\_graph

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph"></a>

## AntiRCGraph Objects

```python
class AntiRCGraph(SchubertMonomialGraph, GridPrint, CrystalGraph)
```

<a id="schubmult.combinatorics.anti_rc_graph.AntiRCGraph.squash_decomp"></a>

#### squash\_decomp

```python
def squash_decomp()
```

Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass.

<a id="schubmult.combinatorics.hpd"></a>

# schubmult.combinatorics.hpd

Bumpless Pipe Dreams (HPD) module

<a id="schubmult.combinatorics.hpd.HPDTile"></a>

## HPDTile Objects

```python
class HPDTile(IntEnum)
```

Enumeration of the possible tile types in a pipe dream.

Each tile represents how two pipes (horizontal and vertical) interact in a square.
Note: Whether a tile is "weighty" depends on the row's _id_vector value, not the tile itself.

<a id="schubmult.combinatorics.hpd.HPDTile.TBD"></a>

#### TBD

Placeholder for uninitialized tile

<a id="schubmult.combinatorics.hpd.HPDTile.BLANK"></a>

#### BLANK

Both pipes go straight (no crossing, no elbow)

<a id="schubmult.combinatorics.hpd.HPDTile.CROSS"></a>

#### CROSS

Pipes cross each other

<a id="schubmult.combinatorics.hpd.HPDTile.HORIZ"></a>

#### HORIZ

Horizontal pipe

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_NW"></a>

#### ELBOW\_NW

Elbow: bottom-right to top-left (╯)

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_SE"></a>

#### ELBOW\_SE

Elbow: top-left to bottom-right (╮)

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_NE"></a>

#### ELBOW\_NE

Elbow: bottom-right to top-left (╯)

<a id="schubmult.combinatorics.hpd.HPDTile.ELBOW_SW"></a>

#### ELBOW\_SW

Elbow: top-left to bottom-right (╮)

<a id="schubmult.combinatorics.hpd.HPDTile.BUMP"></a>

#### BUMP

Bump/osculating tile (pipes touch at corner)

<a id="schubmult.combinatorics.hpd.HPDTile.__str__"></a>

#### \_\_str\_\_

```python
def __str__() -> str
```

Return base display symbol (use get_display_symbol() for context-aware rendering)

<a id="schubmult.combinatorics.hpd.HPDTile.get_display_symbol"></a>

#### get\_display\_symbol

```python
def get_display_symbol(is_weighty: bool) -> str
```

Get display symbol based on whether the tile is in a weighty row

<a id="schubmult.combinatorics.hpd.HPDTile.from_tiletype"></a>

#### from\_tiletype

```python
@classmethod
def from_tiletype(cls, tile: TileType) -> HPDTile
```

Convert from TileType to HPDTile

<a id="schubmult.combinatorics.hpd.HPDTile.is_crossing"></a>

#### is\_crossing

```python
@cached_property
def is_crossing() -> bool
```

True if this tile is a crossing

<a id="schubmult.combinatorics.hpd.HPDTile.is_elbow"></a>

#### is\_elbow

```python
@cached_property
def is_elbow() -> bool
```

True if this tile is any type of elbow

<a id="schubmult.combinatorics.hpd.HPDTile.is_empty"></a>

#### is\_empty

```python
@cached_property
def is_empty() -> bool
```

True if this tile is empty (pipes go straight)

<a id="schubmult.combinatorics.hpd.HPDTile.feeds_right"></a>

#### feeds\_right

```python
@cached_property
def feeds_right() -> bool
```

True if the horizontal pipe continues to the right

<a id="schubmult.combinatorics.hpd.HPDTile.feeds_up"></a>

#### feeds\_up

```python
@cached_property
def feeds_up() -> bool
```

True if the vertical pipe continues upwards

<a id="schubmult.combinatorics.hpd.HPDTile.entrance_from_bottom"></a>

#### entrance\_from\_bottom

```python
@cached_property
def entrance_from_bottom() -> bool
```

True if a pipe can enter from the bottom

<a id="schubmult.combinatorics.hpd.HPDTile.entrance_from_left"></a>

#### entrance\_from\_left

```python
@cached_property
def entrance_from_left() -> bool
```

True if a pipe can enter from the left

<a id="schubmult.combinatorics.hpd.HPD"></a>

## HPD Objects

```python
class HPD(SchubertMonomialGraph, DefaultPrinting)
```

Bumpless Pipe Dream representation.

A bumpless pipe dream is an n×n grid where:
- HPDTile.CROSS (1) represents a crossing
- HPDTile.BLANK (0) represents an empty box (pipes go straight)
- For general pipe dreams, can use HPDTile.ELBOW_* (2-5) for elbows

Each HPD corresponds to a permutation and has an associated weight.

<a id="schubmult.combinatorics.hpd.HPD.__init__"></a>

#### \_\_init\_\_

```python
def __init__(grid, id_vector, *, _is_copy=False) -> None
```

Initialize a HPD from a grid.

**Arguments**:

- `grid` - n×n array-like of HPDTile values, integers 0-5, or list of lists

<a id="schubmult.combinatorics.hpd.HPD.is_classic_row"></a>

#### is\_classic\_row

```python
def is_classic_row(row: int) -> bool
```

Check if row is a classic row (id_vector[row] == 0)

<a id="schubmult.combinatorics.hpd.HPD.concat"></a>

#### concat

```python
@classmethod
def concat(cls, rc, bpd)
```

Concatenate a BPD and RCGraph into a HPD.

**Arguments**:

- `bpd` - BPD instance
- `rc` - RCGraph instance

<a id="schubmult.combinatorics.hpd.HPD.is_weighty_position"></a>

#### is\_weighty\_position

```python
def is_weighty_position(row: int, col: int) -> bool
```

Check if a position is in; a weighty row (id_vector[row] == 1)

<a id="schubmult.combinatorics.hpd.HPD.row_index_to_label"></a>

#### row\_index\_to\_label

```python
def row_index_to_label(row_index: int) -> int
```

Map physical row index (0-based) to row label.

Row labels are assigned counterclockwise:
- _id_vector == 1 rows: labeled 1, 2, ... on RIGHT side, going UPWARD (bottom to top)
- _id_vector == 0 rows: labeled next, on LEFT side, going DOWNWARD (top to bottom)

**Arguments**:

- `row_index` - Physical row index (0-based, top to bottom)
  

**Returns**:

  Row label (1-based)

<a id="schubmult.combinatorics.hpd.HPD.row_label_to_index"></a>

#### row\_label\_to\_index

```python
def row_label_to_index(label: int) -> int
```

Map row label (1-based) to physical row index (0-based).

Inverse of row_index_to_label.

**Arguments**:

- `label` - Row label (1-based)
  

**Returns**:

  Physical row index (0-based, top to bottom)

<a id="schubmult.combinatorics.hpd.HPD.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc: RCGraph) -> HPD
```

Create a HPD from an RC graph.

**Arguments**:

- `rc` - RCGraph instance

<a id="schubmult.combinatorics.hpd.HPD.swap_rows"></a>

#### swap\_rows

```python
def swap_rows(row: int) -> HPD
```

Swap a classic row with the row below it.

**Arguments**:

- `row` - Index of the classic row to swap (0-based)
  

**Returns**:

  New HPD with the specified rows swapped

<a id="schubmult.combinatorics.hpd.HPD.pipe_source_labels"></a>

#### pipe\_source\_labels

```python
def pipe_source_labels(row: int, col: int) -> dict[str, int | None]
```

Determine which row label(s) the pipe(s) at position (row, col) came from.

Returns a dict with keys 'top', 'bottom', 'left', 'right' indicating which
row label the pipe on each side of the tile belongs to (None if no pipe on that side).

Invariants:
- For non-BUMP/CROSS tiles: exactly 2 non-None sides with equal values
- For BUMP: right == bottom, left == top
- For CROSS: top == bottom, left == right

**Arguments**:

- `row` - Physical row index (0-based)
- `col` - Physical column index (0-based)
  

**Returns**:

  Dict with keys 'top', 'bottom', 'left', 'right' mapping to row labels or None

<a id="schubmult.combinatorics.hpd.HPD.from_bruhat_path"></a>

#### from\_bruhat\_path

```python
@classmethod
def from_bruhat_path(cls, path: Sequence[Permutation]) -> HPD
```

Create a HPD from a Bruhat path.

<a id="schubmult.combinatorics.hpd.HPD.row_from_k_chain"></a>

#### row\_from\_k\_chain

```python
@staticmethod
def row_from_k_chain(u: Permutation, w: Permutation, k: int,
                     n: int) -> np.ndarray
```

Construct a single row of tiles from a k-chain according to Definition 3.15.

Given two permutations u and w where u ≤ w in Bruhat order, finds a maximal
k-chain from u to w and constructs a row of n tiles based on that chain.

**Arguments**:

- `u` - Starting permutation
- `w` - Target permutation (must satisfy u ≤ w in Bruhat order)
- `k` - The chain parameter (k >= 1)
  

**Returns**:

  1D numpy array of HPDTile values representing the row
  
  Definition 3.15 cases (for tile at position (row, c)):
  - If chain swaps c with larger but not smaller: ELBOW_SE (⌜)
  - If chain swaps c with both larger and smaller: CROSS (╋)
  - If chain swaps c with smaller but not larger: ELBOW_NW (⌟)
  - If c not among first k numbers of w: BLANK (□)
  - If chain swaps values a,b with a < c < b: BUMP (╬)
  - Otherwise: CROSS (■)

<a id="schubmult.combinatorics.hpd.HPD.build"></a>

#### build

```python
def build() -> None
```

Build internal structures by resolving TBD tiles using lookup table.

<a id="schubmult.combinatorics.hpd.HPD.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Return the size n of the n×n grid

<a id="schubmult.combinatorics.hpd.HPD.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key) -> HPDTile | np.ndarray
```

Access grid elements, casting to HPDTile

<a id="schubmult.combinatorics.hpd.HPD.shiftup"></a>

#### shiftup

```python
def shiftup(shift: int = 1) -> HPD
```

Shift the HPD up by a given amount.

<a id="schubmult.combinatorics.hpd.HPD.perm"></a>

#### perm

```python
@property
def perm() -> Permutation
```

Compute the permutation associated with this HPD.

The permutation is determined by following each vertical pipe from bottom to top.
Pipes enter from the bottom (vertical) and left (horizontal).

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.hpd.HPD.permutation"></a>

#### permutation

```python
@property
def permutation() -> Permutation
```

Alias for perm property

<a id="schubmult.combinatorics.hpd.HPD.inv"></a>

#### inv

```python
@property
def inv() -> int
```

Return the inversion count of the associated permutation.

This is a convenience property that delegates to perm.inv.

**Returns**:

  Number of inversions in the permutation

<a id="schubmult.combinatorics.hpd.HPD.length_vector"></a>

#### length\_vector

```python
@property
def length_vector() -> tuple[int, ...]
```

Compute the length vector of the permutation represented by this HPD.

The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
of weighty tiles in row i. Which tiles are weighty depends on _id_vector[i]:
- If _id_vector[i] == 1: BLANK tiles are weighty
- If _id_vector[i] == 0: CROSS and HORIZ tiles are weighty

**Returns**:

  Tuple of integers representing the length vector

<a id="schubmult.combinatorics.hpd.HPD.from_asm"></a>

#### from\_asm

```python
@classmethod
def from_asm(cls, asm) -> HPD
```

Create a HPD from an ASM (Alternating Sign Matrix).

**Arguments**:

- `asm` - n×n array-like of integers (-1, 0, 1)

**Returns**:

  HPD object

<a id="schubmult.combinatorics.hpd.HPD.weight"></a>

#### weight

```python
@property
def weight() -> Tuple[int, ...]
```

Compute the weight of this HPD.

The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
of empty squares (0s) in column i.

**Returns**:

  Tuple of integers representing the weight

<a id="schubmult.combinatorics.hpd.HPD.word"></a>

#### word

```python
@property
def word() -> Tuple[int, ...]
```

Compute a reduced word for the permutation represented by this HPD.

For each crossing at position (i,j), the word value is the number of pipes
weakly northeast of the crossing minus 1. Weakly northeast means all positions
(r,c) where r <= i and c >= j.

**Returns**:

  Tuple of integers representing the reduced word (1-indexed positions)

<a id="schubmult.combinatorics.hpd.HPD.set_width"></a>

#### set\_width

```python
def set_width(width)
```

Set the width of the HPD by adding empty columns on the right if needed.

<a id="schubmult.combinatorics.hpd.HPD.is_valid"></a>

#### is\_valid

```python
@property
def is_valid() -> bool
```

Check if this is a valid bpd.

**Returns**:

  True if valid, False otherwise

<a id="schubmult.combinatorics.hpd.HPD.__eq__"></a>

#### \_\_eq\_\_

```python
def __eq__(other: object) -> bool
```

Check equality of two HPDs

<a id="schubmult.combinatorics.hpd.HPD.__hash__"></a>

#### \_\_hash\_\_

```python
def __hash__() -> int
```

Hash for use in sets and dicts

<a id="schubmult.combinatorics.hpd.HPD.copy"></a>

#### copy

```python
def copy() -> HPD
```

Create a copy of this HPD

<a id="schubmult.combinatorics.hpd.HPD.num_crossings"></a>

#### num\_crossings

```python
@property
def num_crossings() -> int
```

Total number of crossings in the HPD

<a id="schubmult.combinatorics.hpd.HPD.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the HPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.hpd.HPD.left_root_at"></a>

#### left\_root\_at

```python
def left_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the HPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.hpd.HPD.monk_insert"></a>

#### monk\_insert

```python
def monk_insert(row)
```

RETURNS NORMALIZED

<a id="schubmult.combinatorics.hpd.HPD.combine"></a>

#### combine

```python
def combine(other, shift=None) -> HPD
```

Shift the HPD up by adding empty rows at the bottom.

<a id="schubmult.combinatorics.hpd.HPD.product"></a>

#### product

```python
def product(other: HPD) -> dict[HPD, int]
```

Compute the product of this HPD with another.

<a id="schubmult.combinatorics.hpd.HPD.prod_with_bpd"></a>

#### prod\_with\_bpd

```python
def prod_with_bpd(other: HPD) -> HPD
```

Deprecated: Use product() instead. Returns the single HPD from product dictionary.

<a id="schubmult.combinatorics.hpd.HPD.rebuild"></a>

#### rebuild

```python
def rebuild() -> None
```

Rebuild the HPD to resolve any TBD tiles

<a id="schubmult.combinatorics.hpd.HPD.snap_width"></a>

#### snap\_width

```python
def snap_width() -> HPD
```

Snap the width of the HPD to the length of its permutation.

<a id="schubmult.combinatorics.hpd.HPD.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x: Sequence[Expr],
              y: Sequence[Expr] | None = None,
              **_kwargs) -> Expr
```

Compute the Schubert polynomial value for this HPD.

**Arguments**:

- `x` - Variable or list of variables for polynomial
- `y` - Optional second set of variables for double Schubert polynomial
- ```**_kwargs``` - Additional keyword arguments for polynomial computation (unused)

<a id="schubmult.combinatorics.hpd.HPD.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph() -> RCGraph
```

Convert this HPD to an RC-graph representation.

**Returns**:

  RCGraph object (if available in the module)

<a id="schubmult.combinatorics.schubert_monomial_graph"></a>

# schubmult.combinatorics.schubert\_monomial\_graph

Base class for Schubert monomial graph structures.

Provides a common interface for structures like RCGraph and BPD that represent
monomials in Schubert polynomials via grid-based diagrams.

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph"></a>

## SchubertMonomialGraph Objects

```python
class SchubertMonomialGraph(ABC)
```

Abstract base class for Schubert monomial graph structures.

This class provides a common interface for combinatorial objects that:
- Represent monomials in Schubert polynomials
- Have a grid/matrix-like structure with rows and columns
- Correspond to a permutation
- Have a weight (as a vector or tuple)

Concrete implementations include:
- RCGraph: Reduced-compatible graphs with crystal structure
- BPD: Bumpless pipe dreams

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.perm"></a>

#### perm

```python
@property
@abstractmethod
def perm() -> Permutation
```

Return the permutation associated with this monomial graph.

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.rows"></a>

#### rows

```python
@property
@abstractmethod
def rows() -> int
```

Number of rows in the grid representation.

**Returns**:

  Number of rows

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.cols"></a>

#### cols

```python
@property
@abstractmethod
def cols() -> int
```

Number of columns in the grid representation.

**Returns**:

  Number of columns

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.width"></a>

#### width

```python
@property
def width() -> int
```

Width of the grid (alias for cols).

**Returns**:

  Number of columns

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.height"></a>

#### height

```python
@property
def height() -> int
```

Height of the grid (alias for rows).

**Returns**:

  Number of rows

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.permutation"></a>

#### permutation

```python
@property
def permutation() -> Permutation
```

Alias for perm property.

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.__getitem__"></a>

#### \_\_getitem\_\_

```python
@abstractmethod
def __getitem__(key) -> Any
```

Access elements of the grid.

**Arguments**:

- `key` - Index or tuple of indices
  

**Returns**:

  Element(s) at the specified position

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.normalize"></a>

#### normalize

```python
@abstractmethod
def normalize() -> "SchubertMonomialGraph"
```

Return a normalized version of the monomial graph.

Normalization may involve reordering or simplifying the structure
while preserving its combinatorial properties.

**Returns**:

  Normalized SchubertMonomialGraph object

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.polyvalue"></a>

#### polyvalue

```python
@abstractmethod
def polyvalue(x, y=None, **kwargs) -> Expr
```

Compute the polynomial value represented by this monomial graph.

**Returns**:

  Polynomial representation (e.g., as a SymPy expression)

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.right_zero_act"></a>

#### right\_zero\_act

```python
@abstractmethod
def right_zero_act() -> set["SchubertMonomialGraph"]
```

Compute the right action of a zero (adding a row/column).

**Returns**:

  Set of resulting monomial graphs

<a id="schubmult.combinatorics.schubert_monomial_graph.SchubertMonomialGraph.product"></a>

#### product

```python
@abstractmethod
def product(
        other: "SchubertMonomialGraph") -> dict["SchubertMonomialGraph", int]
```

Compute the product of this monomial graph with another.

This represents the Schubert polynomial multiplication at the monomial level.

**Arguments**:

- `other` - Another monomial graph to multiply with
  

**Returns**:

  Dictionary mapping result monomial graphs to their coefficients

<a id="schubmult.combinatorics.rc_graph"></a>

# schubmult.combinatorics.rc\_graph

<a id="schubmult.combinatorics.rc_graph.RCGraph"></a>

## RCGraph Objects

```python
class RCGraph(WCGraph, CrystalGraph)
```

<a id="schubmult.combinatorics.rc_graph.RCGraph.squash_decomp"></a>

#### squash\_decomp

```python
@cache
def squash_decomp()
```

Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass.

<a id="schubmult.combinatorics.rc_graph.RCGraph.left_squash_decomp"></a>

#### left\_squash\_decomp

```python
def left_squash_decomp()
```

Decompose an n-row RC graph into a pair of n-row RC graph in S_n and an n-grass.

<a id="schubmult.combinatorics.rc_graph.RCGraph.args"></a>

#### args

```python
@property
def args() -> tuple
```

Return args for sympy compatibility - prevents traversal into tuple contents.

<a id="schubmult.combinatorics.rc_graph.RCGraph.extract_demazure_atom"></a>

#### extract\_demazure\_atom

```python
def extract_demazure_atom()
```

Extract the Demazure atom associated with this RC graph.

The Demazure atom is indexed by the right key of the weight tableau,
computed using the Willis algorithm with Earliest Weakly Increasing Subsequences (EWIS).

It consists of all RC graphs with the same permutation and length whose
weight tableaux have the same right key.

Returns a list of RC graphs forming the Demazure atom.

<a id="schubmult.combinatorics.rc_graph.RCGraph.demazure_weight"></a>

#### demazure\_weight

```python
@property
def demazure_weight() -> tuple[int, ...]
```

Weight of the distinguished Demazure extremal element.

<a id="schubmult.combinatorics.rc_graph.RCGraph.w_key_cache"></a>

#### w\_key\_cache

noqa: RUF012

<a id="schubmult.combinatorics.rc_graph.RCGraph.rc_cache"></a>

#### rc\_cache

noqa: RUF012

<a id="schubmult.combinatorics.rc_graph.RCGraph.product"></a>

#### product

```python
@cache
def product(other: RCGraph) -> dict[RCGraph, int]
```

Compute the product of this RC graph with another.

<a id="schubmult.combinatorics.rc_graph.RCGraph.prod_with_rc"></a>

#### prod\_with\_rc

```python
def prod_with_rc(other: RCGraph) -> dict[RCGraph, int]
```

Deprecated: Use product() instead.

<a id="schubmult.combinatorics.rc_graph.RCGraph.demazure_isomorphism_class"></a>

#### demazure\_isomorphism\_class

```python
@property
def demazure_isomorphism_class() -> tuple[tuple[int], Permutation]
```

Alias for classify_demazure_crystal for explicit API usage.

<a id="schubmult.combinatorics.increasing_tableau"></a>

# schubmult.combinatorics.increasing\_tableau

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau"></a>

## IncreasingTableau Objects

```python
class IncreasingTableau(Plactic)
```

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_insert"></a>

#### hecke\_insert

```python
def hecke_insert(*letters)
```

Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.ed_insert_rsk"></a>

#### ed\_insert\_rsk

```python
@classmethod
def ed_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.ed_column_insert_rsk"></a>

#### ed\_column\_insert\_rsk

```python
@classmethod
def ed_column_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this IncreasingTableau tableau and return a new Plactic.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(*corners)
```

K-theoretic (Buch–Samuel) jeu de taquin slide toward the top-left.

Starting from one or more empty inner cells, slide entries up/left into
the holes. This is the genuine K-theoretic slide, so a single entry may
migrate into several holes at once; passing multiple corners performs
the simultaneous slide of all of them.

Corners may be given either as ``up_jdt_slide(row, col)`` (a single
corner) or as ``up_jdt_slide((r1, c1), (r2, c2), ...)``. Returns a new
:class:`IncreasingTableau`.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(*corners)
```

K-theoretic (Buch–Samuel) jeu de taquin slide toward the bottom-right.

Starting from one or more empty outer cells, slide entries down/right
into the holes. As with :meth:`up_jdt_slide`, this is the genuine
K-theoretic slide (an entry may migrate into several holes) and passing
multiple corners performs the simultaneous slide of all of them.

Corners may be given either as ``down_jdt_slide(row, col)`` (a single
corner) or as ``down_jdt_slide((r1, c1), (r2, c2), ...)``. Returns a new
:class:`IncreasingTableau`.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.down_jdt_slide_all_inner_corners"></a>

#### down\_jdt\_slide\_all\_inner\_corners

```python
def down_jdt_slide_all_inner_corners()
```

Simultaneously down-slide every valid inner corner.

Collects all inner corners of the skew shape (see
:attr:`iter_inner_corners`) and performs a single K-theoretic down slide
that moves all of them at once. Returns a new
:class:`IncreasingTableau`. If there are no inner corners the tableau is
returned unchanged.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_column_insert"></a>

#### hecke\_column\_insert

```python
def hecke_column_insert(*letters)
```

Column Hecke-insert ``letters`` into a copy of this tableau.

Returns the resulting :class:`IncreasingTableau` (the insertion tableau
``P``). Use :meth:`hecke_column_insert_rsk` when the set-valued
recording tableau ``Q`` is also required.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_column_insert_rsk"></a>

#### hecke\_column\_insert\_rsk

```python
@classmethod
def hecke_column_insert_rsk(cls, recording, insertion)
```

Column Hecke (K-theoretic) insertion of a two-line array.

The two-line array is given by two equal-length sequences: ``recording``
holds the top-row labels ``k`` (weakly increasing, with the letters in
each equal-label block strictly increasing) and ``insertion`` holds the
bottom-row letters ``a`` that are column-inserted.

Returns ``(P, Q)`` where ``P`` is an :class:`IncreasingTableau` (the
insertion tableau) and ``Q`` is a
:class:`~schubmult.combinatorics.set_valued_tableau.SetValuedTableau`,
the semistandard *set-valued* recording tableau.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.reverse_hecke_column_insert"></a>

#### reverse\_hecke\_column\_insert

```python
def reverse_hecke_column_insert(corner, alpha=1)
```

Reverse Hecke-column-insert the box at ``corner`` ``(row, col)``.

``alpha`` is ``1`` when the corner box was created by the forward
insertion (a ``"grow"`` step) and ``0`` when the forward step merely
absorbed a letter at that corner. Returns ``(Y, x)`` where ``Y`` is the
resulting :class:`IncreasingTableau` and ``x`` the reconstructed letter.

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.hecke_column_uninsert_rsk"></a>

#### hecke\_column\_uninsert\_rsk

```python
@classmethod
def hecke_column_uninsert_rsk(cls, P, Q)
```

Invert :meth:`hecke_column_insert_rsk`.

Given an insertion tableau ``P`` (:class:`IncreasingTableau`) and a
set-valued recording tableau ``Q`` of the same shape, reconstruct the
two-line array as ``(recording, insertion)``.

The recording labels are undone from largest to smallest. Within a box,
the maximum label was the last placed there: a singleton box came from a
``"grow"`` step (reverse with ``alpha = 1``), while a box holding several
labels had its largest label appended by an ``"absorb"`` step (reverse
with ``alpha = 0``, keeping the box). Because insertion is by *columns*
with the canonical convention that equal-label letters are inserted in
strictly *decreasing* order, the corners sharing the maximal label are
undone right-most column first (that being the most recently created
box for the batch).

<a id="schubmult.combinatorics.increasing_tableau.IncreasingTableau.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other)
```

Plactic product: insert entries of `other` in row-reading order
(top-to-bottom, left-to-right) into a copy of self.

<a id="schubmult.combinatorics.set_word"></a>

# schubmult.combinatorics.set\_word

<a id="schubmult.combinatorics.set_word.SetLetter"></a>

## SetLetter Objects

```python
class SetLetter(CrystalGraph, frozenset)
```

A set of ints with a sqrt(gl_n) crystal structure.

<a id="schubmult.combinatorics.set_word.SetWord"></a>

## SetWord Objects

```python
class SetWord(CrystalGraphTensor)
```

A tuple of SetLetters with a sqrt(gl_n) crystal structure.

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

<a id="schubmult.combinatorics.mbpd"></a>

# schubmult.combinatorics.mbpd

Implementation of the MBPD <-> RCP <-> WCGraph bijection from
``writing/mbpd.solve.tex`` (Theorem "T: main").

Geometry (paper convention):
  * n x n grid, rows top->bottom (1..n), columns left->right (1..n).
  * Pipes ENTER from the RIGHT border and EXIT to the SOUTH border.
  * No pipe enters from the top, no pipe leaves to the left.

Seven tiles, described by the set of directions in which the pipe connects
(N=up, E=right, S=down, W=left) together with a "marked" bit:

    B  blank        {}            heavy
    H  horizontal   {E, W}
    V  vertical     {N, S}
    P  plus/cross   {N, E, S, W}
    R  R-elbow      {E, S}
    J  J-elbow      {N, W}
    M  marked J     {N, W}        heavy

A tile is *heavy* iff it is B or M.

The bijection ``Phi = MBPD.phi`` and its inverse ``Psi = RCP.psi`` preserve the
associated permutation and weight.  ``WCGraph.to_mbpd`` / ``WCGraph.from_mbpd``
expose the composition with the trivial ``RCP <-> WCGraph`` repackaging.

<a id="schubmult.combinatorics.mbpd.tile_name"></a>

#### tile\_name

```python
def tile_name(conn: frozenset, marked: bool) -> str
```

Return the canonical tile name for a connection set / mark.

<a id="schubmult.combinatorics.mbpd.MBPD"></a>

## MBPD Objects

```python
class MBPD()
```

A marked bumpless pipedream on an ``n x n`` grid.

Internally we store, for every cell ``(i, j)`` (1-indexed), the connection
set ``conn[i][j]`` (a frozenset of ``N/E/S/W``) and a boolean ``marked``.

<a id="schubmult.combinatorics.mbpd.MBPD.from_tiles"></a>

#### from\_tiles

```python
@classmethod
def from_tiles(cls, grid) -> MBPD
```

Build from an ``n x n`` grid of tile-name strings.

<a id="schubmult.combinatorics.mbpd.MBPD.with_tile"></a>

#### with\_tile

```python
def with_tile(i: int, j: int, name: str) -> MBPD
```

Return a copy with cell ``(i, j)`` set to tile ``name``.

<a id="schubmult.combinatorics.mbpd.MBPD.rothe"></a>

#### rothe

```python
@classmethod
def rothe(cls, perm, n=None) -> MBPD
```

The Rothe MBPD ``D_w``: pipe ``i`` turns only at ``(i, w(i))``.

Connections (derived and verified against the paper's conventions):
  E : j >= w(i)                 (horizontal present / border on right)
  W : j > w(i)
  S : w^{-1}(j) <= i            (vertical present going down)
  N : w^{-1}(j) < i

<a id="schubmult.combinatorics.mbpd.MBPD.validity_errors"></a>

#### validity\_errors

```python
def validity_errors()
```

Return a list of human-readable validity problems (empty if valid).

<a id="schubmult.combinatorics.mbpd.MBPD.perm"></a>

#### perm

```python
def perm()
```

Associated permutation ``w`` where the pipe entering the right of
row ``i`` leaves the bottom of column ``w(i)``.

Double crossings between an already-crossed pair are ignored (the
Demazure convention): the *first* time a pair of pipes meets at a ``P``
tile they cross, every later meeting is a *bump* (the ``P`` acts as
superimposed elbows ``R``+``J``, i.e. ``E<->S`` and ``N<->W``).

We compute the routing by a fixed point.  Start assuming every ``P`` is
a genuine crossing (straight through) and trace all pipes.  Whenever a
pair of pipes meets at more than one ``P`` tile we mark all but their
first meeting as bumps, then re-trace.  This is repeated until the set
of bump tiles stabilises; the exit columns then give ``w``.

<a id="schubmult.combinatorics.mbpd.MBPD.weight"></a>

#### weight

```python
def weight()
```

``wt(D) = (m_1, ..., m_n)`` with ``m_i`` = `heavy` tiles in row ``i``.

<a id="schubmult.combinatorics.mbpd.MBPD.is_pipe_segment"></a>

#### is\_pipe\_segment

```python
def is_pipe_segment(r: int, b: int, c: int) -> bool
```

``D_{r,[b,c]}`` is a pipe segment: a connected horizontal run.

For ``b == c`` a single non-blank tile qualifies.  For ``b < c`` the
interior ``D_{r,[b+1,c-1]}`` must be all ``H``/``P`` tiles and the run
must actually connect: ``(r,b)`` connects east and ``(r,c)`` connects
west (so the paper's inference that the endpoints connect toward the
interior holds).

<a id="schubmult.combinatorics.mbpd.MBPD.rj_subsequence"></a>

#### rj\_subsequence

```python
def rj_subsequence(r: int, lo: int, hi: int)
```

The RJ subsequence (list of 'R'/'J') of light tiles in ``[lo,hi]``.

Returns ``None`` if any tile in the range is heavy (not a light seq).

<a id="schubmult.combinatorics.mbpd.MBPD.light_seq_type"></a>

#### light\_seq\_type

```python
def light_seq_type(r: int, lo: int, hi: int)
```

Classify the light sequence ``D_{r,[lo,hi]}`` as one of
``'paired'``, ``'J'``, ``'R'``, ``'JR'`` (or ``None`` if it is not a
light sequence).  Empty range counts as ``'paired'``.

<a id="schubmult.combinatorics.mbpd.MBPD.is_doublecross"></a>

#### is\_doublecross

```python
def is_doublecross(r: int, b: int, d: int) -> bool
```

``D_{[r,r+1],[b,d]}`` is a doublecross: both rows are pipe segments,
``D_{r,b}=R`` and ``D_{r+1,d}=J``.

<a id="schubmult.combinatorics.mbpd.MBPD.max_F_target"></a>

#### max\_F\_target

```python
def max_F_target()
```

Bottommost then rightmost heavy tile, or ``None`` for ``D_id``.

<a id="schubmult.combinatorics.mbpd.MBPD.f_move"></a>

#### f\_move

```python
def f_move(r: int, c: int) -> MBPD
```

Apply the ``F``-move at the ``F``-target ``(r,c)`` (paper
\S "do F move").

<a id="schubmult.combinatorics.mbpd.MBPD.e_target"></a>

#### e\_target

```python
def e_target(r: int)
```

Determine the ``E``-target ``(r+1,c)`` in the two-row window with top
row ``r`` (``1 <= r < n``).  Returns a dict with keys
``kind`` (``'e'`` or ``'estar'``), ``c``, ``cprime`` or ``None`` if no
``E``-target exists for this row.

<a id="schubmult.combinatorics.mbpd.MBPD.E_target_info"></a>

#### E\_target\_info

```python
def E_target_info(r: int) -> dict
```

Full case analysis for the ``E``-move at row ``r``: the right
trichotomy (Initial/Plus/NPlus) with ``rho``, and the left trichotomy
(Straight/DCross/Leftturn) with ``lambda``.

<a id="schubmult.combinatorics.mbpd.MBPD.e_move"></a>

#### e\_move

```python
def e_move(r: int) -> MBPD
```

Apply the ``E``-move ``E_r`` at row ``r`` (paper \S "the E moves"):
* if ``(r+1,c)`` is ``M``, unmark it;
* in Cases Straight and DCross do the ``(r,[lambda,rho])``-droop;
* if ``(r,lambda)`` is ``J``, mark it.

<a id="schubmult.combinatorics.mbpd.MBPD.phi"></a>

#### phi

```python
def phi() -> RCP
```

The bijection ``Phi`` of the paper (Theorem "T: main"), realized by
the row-pop recursion:

  * ``Phi(D_id) = ()``
  * F-terminal:    ``Phi(D) = (i,i) . Phi(f*_i(D))``
  * F-nonterminal: ``Phi(D) = down( Phi(f_i(D)) )``

where ``i`` is the row of the maximum ``F``-target and ``down`` replaces
the first biletter ``(i',a)`` by ``(i'-1,a)`` (RRCP2).

<a id="schubmult.combinatorics.mbpd.MBPD.to_bpd"></a>

#### to\_bpd

```python
def to_bpd()
```

Return ``(bpd, marks)``: the underlying :class:`BPD` and the frozenset
of ``(i, j)`` cells (1-indexed) that carry a mark (``M`` tiles).

<a id="schubmult.combinatorics.mbpd.MBPD.from_bpd"></a>

#### from\_bpd

```python
@classmethod
def from_bpd(cls, bpd, marks=frozenset()) -> MBPD
```

Rebuild an :class:`MBPD` from a :class:`BPD` and a set of marked
``(i, j)`` cells.  The BPD is squared up to ``max(rows, len(perm))`` first
(Rothe completion), and a mark is only reinstated where the tile is a
``J`` up-elbow.

<a id="schubmult.combinatorics.mbpd.MBPD.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
def zero_out_last_row(rows: int) -> MBPD
```

Zero out below ``rows`` by resizing the underlying BPD down to ``rows``
rows (Weigandt/Lascoux transition) and re-deducing the marks that survive.

Weight is preserved; the permutation changes according to the transition.

<a id="schubmult.combinatorics.mbpd.RCP"></a>

## RCP Objects

```python
@dataclass(frozen=True)
class RCP()
```

A reverse compatible pair: a tuple of biletters ``(i, a)`` with
``1 <= i <= a < n`` (``n`` is the ambient size), strictly *decreasing* in
the order ``(i1,a1) > (i2,a2)`` iff ``i1 > i2`` or (``i1 == i2`` and
``a1 < a2``).

<a id="schubmult.combinatorics.mbpd.RCP.biletters"></a>

#### biletters

tuple of (i, a)

<a id="schubmult.combinatorics.mbpd.RCP.perm"></a>

#### perm

```python
def perm()
```

``w(B) = s_{a_l} * ... * s_{a_1}`` (Demazure product, subscripts
decreasing = reversed biletter order).

<a id="schubmult.combinatorics.mbpd.RCP.to_wcgraph"></a>

#### to\_wcgraph

```python
def to_wcgraph()
```

Package the RCP as a :class:`WCGraph`.

``WCGraph.from_word_compatible`` wants a *weakly increasing* compatible
sequence.  Reversing the biletter list turns the RCP's decreasing order
into weakly increasing ``i`` with strictly decreasing ``a`` within a
block -- exactly the WCGraph compatibility condition.

<a id="schubmult.combinatorics.mbpd.RCP.from_wcgraph"></a>

#### from\_wcgraph

```python
@classmethod
def from_wcgraph(cls, wc, n=None)
```

Inverse of :meth:`to_wcgraph`.

<a id="schubmult.combinatorics.mbpd.RCP.to_pd_crossings"></a>

#### to\_pd\_crossings

```python
def to_pd_crossings()
```

Remark "R:CP and PD": biletter ``(i, a)`` -> crossing at
``(i, a - i + 1)`` in the ordinary pipedream picture.

<a id="schubmult.combinatorics.mbpd.RCP.psi"></a>

#### psi

```python
def psi() -> MBPD
```

The inverse bijection ``Psi = Phi^{-1}`` (Theorem "T: main"),
realized by row-unpop.

``Phi(D) = (i,a) . Phi(nabla_r D)`` where the maximum ``F``-target of
``D`` is in row ``i`` and ``nabla_r D = f*_a f_{a-1} ... f_i (D)`` is
obtained by climbing rows ``i, i+1, ..., a`` (the last, at row ``a``,
being the ``f*``-move).  Inverting: given the first biletter ``(i,a)``,
first rebuild ``N = nabla_r D = Psi(rest)``, then apply
``D = e_i(e_{i+1}( ... e_{a-1}( e*_a(N) ) ... ))``.

<a id="schubmult.combinatorics.crystal_graph"></a>

# schubmult.combinatorics.crystal\_graph

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph"></a>

## CrystalGraph Objects

```python
class CrystalGraph(Printable)
```

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(index)
```

The raising operator for the crystal graph.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(index)
```

The lowering operator for the crystal graph.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

The weight of the crystal graph element.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.to_lowest_weight"></a>

#### to\_lowest\_weight

```python
def to_lowest_weight(length=None)
```

Return the lowest weight element in the connected component.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.to_highest_weight"></a>

#### to\_highest\_weight

```python
def to_highest_weight(length=None, start=1)
```

Return the highest weight element in the connected component.

<a id="schubmult.combinatorics.crystal_graph.CrystalGraph.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Return the length of the crystal element.

<a id="schubmult.combinatorics.bpd"></a>

# schubmult.combinatorics.bpd

Bumpless Pipe Dreams (BPD) module

<a id="schubmult.combinatorics.bpd.TileType"></a>

## TileType Objects

```python
class TileType(IntEnum)
```

Enumeration of the 6 possible tile types in a pipe dream.

Each tile represents how two pipes (horizontal and vertical) interact in a square.

<a id="schubmult.combinatorics.bpd.TileType.TBD"></a>

#### TBD

Placeholder for uninitialized tile

<a id="schubmult.combinatorics.bpd.TileType.BLANK"></a>

#### BLANK

Both pipes go straight (no crossing, no elbow)

<a id="schubmult.combinatorics.bpd.TileType.CROSS"></a>

#### CROSS

Pipes cross each other

<a id="schubmult.combinatorics.bpd.TileType.ELBOW_NW"></a>

#### ELBOW\_NW

Elbow: bottom-right to top-left (╯)

<a id="schubmult.combinatorics.bpd.TileType.ELBOW_SE"></a>

#### ELBOW\_SE

Elbow: top-left to bottom-right (╮)

<a id="schubmult.combinatorics.bpd.TileType.BUMP"></a>

#### BUMP

Bump/osculating tile (pipes touch at corner)

<a id="schubmult.combinatorics.bpd.TileType.as_tile"></a>

#### as\_tile

```python
def as_tile() -> Tile
```

Convert TileType to a Tile object with edge information.

<a id="schubmult.combinatorics.bpd.TileType.is_crossing"></a>

#### is\_crossing

```python
@cached_property
def is_crossing() -> bool
```

True if this tile is a crossing

<a id="schubmult.combinatorics.bpd.TileType.is_elbow"></a>

#### is\_elbow

```python
@cached_property
def is_elbow() -> bool
```

True if this tile is any type of elbow

<a id="schubmult.combinatorics.bpd.TileType.is_empty"></a>

#### is\_empty

```python
@cached_property
def is_empty() -> bool
```

True if this tile is empty (pipes go straight)

<a id="schubmult.combinatorics.bpd.TileType.feeds_right"></a>

#### feeds\_right

```python
@cached_property
def feeds_right() -> bool
```

True if the horizontal pipe continues to the right

<a id="schubmult.combinatorics.bpd.TileType.feeds_up"></a>

#### feeds\_up

```python
@cached_property
def feeds_up() -> bool
```

True if the vertical pipe continues upwards

<a id="schubmult.combinatorics.bpd.TileType.entrance_from_bottom"></a>

#### entrance\_from\_bottom

```python
@cached_property
def entrance_from_bottom() -> bool
```

True if a pipe can enter from the bottom

<a id="schubmult.combinatorics.bpd.TileType.entrance_from_left"></a>

#### entrance\_from\_left

```python
@cached_property
def entrance_from_left() -> bool
```

True if a pipe can enter from the left

<a id="schubmult.combinatorics.bpd.BPD"></a>

## BPD Objects

```python
class BPD(SchubertMonomialGraph, DefaultPrinting)
```

Bumpless Pipe Dream representation.

A bumpless pipe dream is an n×n grid where:
- TileType.CROSS (1) represents a crossing
- TileType.BLANK (0) represents an empty box (pipes go straight)
- For general pipe dreams, can use TileType.ELBOW_* (2-5) for elbows

Each BPD corresponds to a permutation and has an associated weight.

<a id="schubmult.combinatorics.bpd.BPD.__init__"></a>

#### \_\_init\_\_

```python
def __init__(grid,
             column_perm: Permutation | None = None,
             *,
             _is_copy=False) -> None
```

Initialize a BPD from a grid.

**Arguments**:

- `grid` - n×n array-like of TileType values, integers 0-5, or list of lists

<a id="schubmult.combinatorics.bpd.BPD.all_unreduced_bpds"></a>

#### all\_unreduced\_bpds

```python
@classmethod
def all_unreduced_bpds(cls,
                       w: Permutation,
                       length: int | None = None,
                       weight: tuple[int] | None = None) -> set[BPD]
```

Enumerate all (possibly unreduced) BPDs for w.

Defers to :meth:`WCGraph.all_wc_graphs`, converts each WCGraph to its
marked bumpless pipe dream, forgets the marks, and reduces to an
ordinary BPD. This is dramatically faster than enumerating ASMs.

<a id="schubmult.combinatorics.bpd.BPD.from_bruhat_path"></a>

#### from\_bruhat\_path

```python
@classmethod
def from_bruhat_path(cls, path: Sequence[Permutation]) -> BPD
```

Create a BPD from a Bruhat path.

<a id="schubmult.combinatorics.bpd.BPD.row_from_k_chain"></a>

#### row\_from\_k\_chain

```python
@staticmethod
def row_from_k_chain(u: Permutation, w: Permutation, k: int,
                     n: int) -> np.ndarray
```

Construct a single row of tiles from a k-chain according to Definition 3.15.

Given two permutations u and w where u ≤ w in Bruhat order, finds a maximal
k-chain from u to w and constructs a row of n tiles based on that chain.

**Arguments**:

- `u` - Starting permutation
- `w` - Target permutation (must satisfy u ≤ w in Bruhat order)
- `k` - The chain parameter (k >= 1)
  

**Returns**:

  1D numpy array of TileType values representing the row
  
  Definition 3.15 cases (for tile at position (row, c)):
  - If chain swaps c with larger but not smaller: ELBOW_SE (⌜)
  - If chain swaps c with both larger and smaller: CROSS (╋)
  - If chain swaps c with smaller but not larger: ELBOW_NW (⌟)
  - If c not among first k numbers of w: BLANK (□)
  - If chain swaps values a,b with a < c < b: BUMP (╬)
  - Otherwise: CROSS (■)

<a id="schubmult.combinatorics.bpd.BPD.build"></a>

#### build

```python
def build() -> None
```

Build internal structures by resolving TBD tiles using lookup table.

<a id="schubmult.combinatorics.bpd.BPD.__len__"></a>

#### \_\_len\_\_

```python
def __len__() -> int
```

Return the size n of the n×n grid

<a id="schubmult.combinatorics.bpd.BPD.__getitem__"></a>

#### \_\_getitem\_\_

```python
def __getitem__(key) -> TileType | np.ndarray
```

Access grid elements, casting to TileType

<a id="schubmult.combinatorics.bpd.BPD.shiftup"></a>

#### shiftup

```python
def shiftup(shift: int = 1) -> BPD
```

Shift the BPD up by a given amount.

<a id="schubmult.combinatorics.bpd.BPD.perm"></a>

#### perm

```python
@property
def perm() -> Permutation
```

Compute the permutation associated with this BPD.

The permutation is determined by following each vertical pipe from bottom to top.
Pipes enter from the bottom (vertical) and left (horizontal).

**Returns**:

  Permutation object

<a id="schubmult.combinatorics.bpd.BPD.permutation"></a>

#### permutation

```python
@property
def permutation() -> Permutation
```

Alias for perm property

<a id="schubmult.combinatorics.bpd.BPD.disjoint_union"></a>

#### disjoint\_union

```python
def disjoint_union(other: BPD) -> BPD
```

Row-preserving disjoint union.

This keeps both summands on the same row indices by placing them side-by-side,
with a single horizontal connector column between them.
The result is intended to preserve row placement of blanks and may be unreduced.

<a id="schubmult.combinatorics.bpd.BPD.disjoint_union_block_diag"></a>

#### disjoint\_union\_block\_diag

```python
def disjoint_union_block_diag(other: BPD) -> BPD
```

Disjoint union via ASM block diagonal (rows of the second summand are shifted).

<a id="schubmult.combinatorics.bpd.BPD.disjoint_union_many"></a>

#### disjoint\_union\_many

```python
@classmethod
def disjoint_union_many(cls, *bpds: BPD) -> BPD
```

Row-preserving disjoint union of any number of BPDs.

<a id="schubmult.combinatorics.bpd.BPD.inv"></a>

#### inv

```python
@property
def inv() -> int
```

Return the inversion count of the associated permutation.

This is a convenience property that delegates to perm.inv.

**Returns**:

  Number of inversions in the permutation

<a id="schubmult.combinatorics.bpd.BPD.length_vector"></a>

#### length\_vector

```python
@property
def length_vector() -> tuple[int, ...]
```

Compute the length vector of the permutation represented by this BPD.

The length vector is a tuple (l_1, l_2, ..., l_n) where l_i is the number
of crossings in row i.

**Returns**:

  Tuple of integers representing the length vector

<a id="schubmult.combinatorics.bpd.BPD.from_asm"></a>

#### from\_asm

```python
@classmethod
def from_asm(cls, asm) -> BPD
```

Create a BPD from an ASM (Alternating Sign Matrix).

**Arguments**:

- `asm` - n×n array-like of integers (-1, 0, 1)

**Returns**:

  BPD object

<a id="schubmult.combinatorics.bpd.BPD.weight"></a>

#### weight

```python
@property
def weight() -> Tuple[int, ...]
```

Compute the weight of this BPD.

The weight is a tuple (w_1, w_2, ..., w_n) where w_i is the number
of empty squares (0s) in column i.

**Returns**:

  Tuple of integers representing the weight

<a id="schubmult.combinatorics.bpd.BPD.word"></a>

#### word

```python
@property
def word() -> Tuple[int, ...]
```

Compute a reduced word for the permutation represented by this BPD.

For each crossing at position (i,j), the word value is the number of pipes
weakly northeast of the crossing minus 1. Weakly northeast means all positions
(r,c) where r <= i and c >= j.

**Returns**:

  Tuple of integers representing the reduced word (1-indexed positions)

<a id="schubmult.combinatorics.bpd.BPD.set_width"></a>

#### set\_width

```python
def set_width(width)
```

Set the width of the BPD by adding empty columns on the right if needed.

<a id="schubmult.combinatorics.bpd.BPD.is_valid"></a>

#### is\_valid

```python
@property
def is_valid() -> bool
```

Check if this is a valid bpd.

**Returns**:

  True if valid, False otherwise

<a id="schubmult.combinatorics.bpd.BPD.__eq__"></a>

#### \_\_eq\_\_

```python
def __eq__(other: object) -> bool
```

Check equality of two BPDs

<a id="schubmult.combinatorics.bpd.BPD.__hash__"></a>

#### \_\_hash\_\_

```python
def __hash__() -> int
```

Hash for use in sets and dicts

<a id="schubmult.combinatorics.bpd.BPD.copy"></a>

#### copy

```python
def copy() -> BPD
```

Create a copy of this BPD

<a id="schubmult.combinatorics.bpd.BPD.num_crossings"></a>

#### num\_crossings

```python
@property
def num_crossings() -> int
```

Total number of crossings in the BPD

<a id="schubmult.combinatorics.bpd.BPD.right_root_at"></a>

#### right\_root\_at

```python
def right_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the BPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.bpd.BPD.left_root_at"></a>

#### left\_root\_at

```python
def left_root_at(i: int, j: int) -> int
```

Compute the inversion associated with the crossing at position (i, j).

The inversion is determined by tracing the pipes through the BPD.

**Arguments**:

- `i` - Row index of the crossing
- `j` - Column index of the crossing

**Returns**:

  The inversion value as an integer

<a id="schubmult.combinatorics.bpd.BPD.monk_insert"></a>

#### monk\_insert

```python
def monk_insert(row)
```

RETURNS NORMALIZED

<a id="schubmult.combinatorics.bpd.BPD.combine"></a>

#### combine

```python
def combine(other, shift=None) -> BPD
```

Shift the BPD up by adding empty rows at the bottom.

<a id="schubmult.combinatorics.bpd.BPD.product"></a>

#### product

```python
def product(other: BPD) -> dict[BPD, int]
```

Compute the product of this BPD with another.

<a id="schubmult.combinatorics.bpd.BPD.prod_with_bpd"></a>

#### prod\_with\_bpd

```python
def prod_with_bpd(other: BPD) -> BPD
```

Deprecated: Use product() instead. Returns the single BPD from product dictionary.

<a id="schubmult.combinatorics.bpd.BPD.rebuild"></a>

#### rebuild

```python
def rebuild() -> None
```

Rebuild the BPD to resolve any TBD tiles

<a id="schubmult.combinatorics.bpd.BPD.snap_width"></a>

#### snap\_width

```python
def snap_width() -> BPD
```

Snap the width of the BPD to the length of its permutation.

<a id="schubmult.combinatorics.bpd.BPD.polyvalue"></a>

#### polyvalue

```python
def polyvalue(x: Sequence[Expr],
              y: Sequence[Expr] | None = None,
              **_kwargs) -> Expr
```

Compute the Schubert polynomial value for this BPD.

**Arguments**:

- `x` - Variable or list of variables for polynomial
- `y` - Optional second set of variables for double Schubert polynomial
- ```**_kwargs``` - Additional keyword arguments for polynomial computation (unused)

<a id="schubmult.combinatorics.bpd.BPD.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph() -> RCGraph
```

Convert this BPD to an RC-graph representation.

**Returns**:

  RCGraph object (if available in the module)

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

<a id="schubmult.combinatorics.nilplactic"></a>

# schubmult.combinatorics.nilplactic

<a id="schubmult.combinatorics.nilplactic.NilPlactic"></a>

## NilPlactic Objects

```python
class NilPlactic(Plactic)
```

<a id="schubmult.combinatorics.nilplactic.NilPlactic.all_skew_ed_tableaux"></a>

#### all\_skew\_ed\_tableaux

```python
@classmethod
@cache
def all_skew_ed_tableaux(cls, outer_shape, bruhat_perm, inner_shape=None)
```

Return all skew Edelman–Greene tableaux of skew shape outer_shape/inner_shape
(represented as row lengths) that are increasing (strictly increasing rows
left-to-right and columns top-to-bottom) and whose row-reading word (E-G
row word: read rows bottom-to-top, left-to-right, omitting zeros) has
Permutation.ref_product(...) less than or equal to ``bruhat_perm`` in
Bruhat order.

Representation details:

- outer_shape: sequence of row lengths (one int per row).
- inner_shape: optional sequence of left offsets per row (defaults to zeros).
- The returned tableaux are NilPlactic instances whose rows are left-aligned
  of length ``outer_shape[r]``, with the first ``inner_shape[r]`` entries set
  to 0 to represent the inner (removed) cells of the skew shape.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(row, col)
```

Perform an upward K-theoretic jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new Plactic tableau.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(row, col)
```

Perform a K-theoretic jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new NilPlactic tableau.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.ed_insert"></a>

#### ed\_insert

```python
def ed_insert(*letters)
```

Insert a letter/entry into this NilPlactic tableau and return a new Plactic.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.ed_insert_rsk"></a>

#### ed\_insert\_rsk

```python
@classmethod
def ed_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this NilPlactic tableau and return a new Plactic.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.ed_column_insert_rsk"></a>

#### ed\_column\_insert\_rsk

```python
@classmethod
def ed_column_insert_rsk(cls, w1, w2)
```

Insert a letter/entry into this NilPlactic tableau and return a new Plactic.

<a id="schubmult.combinatorics.nilplactic.NilPlactic.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other)
```

Plactic product: insert entries of `other` in row-reading order
(top-to-bottom, left-to-right) into a copy of self.

<a id="schubmult.combinatorics.root_tableau"></a>

# schubmult.combinatorics.root\_tableau

<a id="schubmult.combinatorics.root_tableau.RootTableau"></a>

## RootTableau Objects

```python
class RootTableau(CrystalGraph, GridPrint)
```

Root tableau with dual knuth equivalence

<a id="schubmult.combinatorics.root_tableau.RootTableau.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(row, col, check=False)
```

Perform a downward/rightward jeu-de-taquin slide starting from the given
(row, col) hole (0-indexed). Boxes from below or to the right are moved
into the hole, preferring the smaller recording letter when both exist.
Returns a new RootTableau (does not mutate self).

<a id="schubmult.combinatorics.root_tableau.RootTableau.rc_graph"></a>

#### rc\_graph

```python
@property
def rc_graph()
```

Reconstruct the RC-graph from the root tableau.

<a id="schubmult.combinatorics.root_tableau.RootTableau.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

Crystal raising operator e_i on the root tableau

<a id="schubmult.combinatorics.root_tableau.RootTableau.raising_operator_direct"></a>

#### raising\_operator\_direct

```python
def raising_operator_direct(i)
```

Direct crystal raising operator e_i using EG invariant tracking.

Key insight: EG invariant is preserved, but eg_index_word changes.

<a id="schubmult.combinatorics.root_tableau.RootTableau.lowering_operator_direct"></a>

#### lowering\_operator\_direct

```python
def lowering_operator_direct(i)
```

Direct crystal lowering operator f_i using EG invariant tracking.

<a id="schubmult.combinatorics.quasi_crystal_graph"></a>

# schubmult.combinatorics.quasi\_crystal\_graph

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph"></a>

## QuasiCrystalGraph Objects

```python
class QuasiCrystalGraph(CrystalGraph)
```

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.quasi_raising_operator"></a>

#### quasi\_raising\_operator

```python
def quasi_raising_operator(index)
```

The raising operator for the crystal graph.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.quasi_lowering_operator"></a>

#### quasi\_lowering\_operator

```python
def quasi_lowering_operator(index)
```

The lowering operator for the crystal graph.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.to_quasi_lowest_weight"></a>

#### to\_quasi\_lowest\_weight

```python
def to_quasi_lowest_weight()
```

Return the lowest weight element in the same quasi-crystal component.

<a id="schubmult.combinatorics.quasi_crystal_graph.QuasiCrystalGraph.reverse_quasi_lower_seq"></a>

#### reverse\_quasi\_lower\_seq

```python
def reverse_quasi_lower_seq(seq)
```

Apply the reverse of the given lowering sequence.

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

<a id="schubmult.combinatorics.wc_graph"></a>

# schubmult.combinatorics.wc\_graph

<a id="schubmult.combinatorics.wc_graph.WCGraph"></a>

## WCGraph Objects

```python
class WCGraph(SchubertMonomialGraph, CrystalGraph, GridPrint, tuple)
```

Word-compatible graph.

Internal representation matches RCGraph: a tuple of rows where row i (0-indexed)
is a strictly decreasing tuple of positive integers >= i + 1.

For WCGraph, the associated permutation is the Demazure product of the
concatenated row word (1-indexed simple reflections).

<a id="schubmult.combinatorics.wc_graph.WCGraph.to_mbpd"></a>

#### to\_mbpd

```python
@cache
def to_mbpd(n: int | None = None)
```

The marked bumpless pipedream ``Psi(RCP(self))`` (paper
``writing/mbpd.solve.tex``, Theorem "T: main").

This composes the trivial ``WCGraph -> RCP`` repackaging with the
row-unpop bijection ``Psi``.  ``n`` is the ambient grid size (defaults
to ``len(self.perm) - 1`` padded to fit the graph).

Cached: WCGraphs are immutable and hashable, so the (self, n) round
trip is memoized.

<a id="schubmult.combinatorics.wc_graph.WCGraph.from_mbpd"></a>

#### from\_mbpd

```python
@classmethod
@cache
def from_mbpd(cls, mbpd) -> WCGraph
```

Inverse of :meth:`to_mbpd`: the WCGraph ``RCP(Phi(mbpd))`` obtained
from the row-pop bijection ``Phi`` (paper Theorem "T: main").

Cached on the (hashable) ``mbpd`` argument.

<a id="schubmult.combinatorics.wc_graph.WCGraph.grove_wcs"></a>

#### grove\_wcs

```python
@classmethod
def grove_wcs(cls, comp, length=None, base_rc=None) -> set[WCGraph]
```

Grove polynomial of ``comp`` built by inverting the omega insertion.

We enumerate the compatible set-valued labelings ``kappa`` of the indexed
forest ``F = weak_composition_to_indfor(comp)`` (the grove definition), and
realize each labeling as a WCGraph by writing down its (word, compatible
sequence) pair *explicitly* -- the published inverse of the set-valued omega
insertion -- and feeding it to ``WCGraph.from_word_compatible``.

Concretely:

* The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
  labeling ``Q`` of the principal reduced RC graph form the omega pair for
  ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
  reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
  so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
  word position ``len(W) - Q(v)`` and carries the reduced letter
  ``ell(v) = W[len(W) - Q(v)]``.
* A labeling ``kappa`` places the letter ``ell(v)`` into every row
  ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
  ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
  a word whose rows are strictly decreasing, i.e. exactly the data
  ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
  values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
  equals ``x^kappa``.

Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
the degree ``-1`` homogenizer recording extra labels beyond one per node.

<a id="schubmult.combinatorics.wc_graph.WCGraph.grove_invariant"></a>

#### grove\_invariant

```python
@property
def grove_invariant() -> tuple[int, ...]
```

The composition of the grove containing ``self`` -- the inverse of
:meth:`grove_wcs`.

Every WCGraph of a fixed permutation lies in exactly one grove.  The grove
is generated by the set-valued labelings of a single indexed forest ``F``,
and its *base* member is the distinguished graph whose ``forest_weight``
equals its ``length_vector`` (the minimal, one-element-per-node labeling).
The grove weight is that base's ``forest_weight`` -- i.e. ``F.code`` padded
to ``len(self)``.

For a base member (``forest_weight == length_vector``) the grove weight
coincides with the forest weight; for set-valued members it recovers the
composition of the base rather than the (larger, occurrence-inflated)
forest weight of ``self`` itself.

<a id="schubmult.combinatorics.wc_graph.WCGraph.strong_hecke_invariant"></a>

#### strong\_hecke\_invariant

```python
@property
def strong_hecke_invariant()
```

K-theoretic rectification of the diagonal-strip increasing tableau.

Lay the ``perm_word`` out along a single anti-diagonal strip, reading
from bottom-left to top-right, as a skew
:class:`~schubmult.combinatorics.increasing_tableau.IncreasingTableau`,
then repeatedly apply the simultaneous inner-corner K-theoretic down
slide until the tableau is rectified (no inner corners remain). The
rectified increasing tableau is the strong Hecke invariant.

<a id="schubmult.combinatorics.wc_graph.WCGraph.to_rc_pieri"></a>

#### to\_rc\_pieri

```python
def to_rc_pieri()
```

Map this WCGraph to an RCGraph via ``_snap_reduced`` + Pieri insertion.

Snap to the underlying reduced RCGraph, then reinsert the missing
(excess) letters row-by-row: the number of letters missing from row
``i + 1`` is ``self.length_vector[i] - reduced.length_vector[i]``, and
those rows (with multiplicity) are Pieri-inserted at ``perm.max_descent``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.zero_out_last_row"></a>

#### zero\_out\_last\_row

```python
@cache
def zero_out_last_row() -> WCGraph
```

Zero out the (empty) last row, realizing the Weigandt/Lascoux
transition (``writing/weigandt_bumpless.tex``, Theorem "transition").

The graph is transported to a marked bumpless pipedream, the underlying
BPD is resized down by one row (dropping the maximal-corner row), and the
result is transported back.  Weight is preserved; the permutation changes
according to the transition.  Unlike the previous Hecke-insertion route,
this is total: it works for non-core-reduced graphs as well.

Cached: the full MBPD round trip is memoized on the (hashable) graph,
so repeated ``squash_product`` calls reuse the result.

<a id="schubmult.combinatorics.wc_graph.WCGraph.product"></a>

#### product

```python
def product(other: SchubertMonomialGraph) -> dict[WCGraph, int]
```

Compute the product of this WC graph with another.

<a id="schubmult.combinatorics.wc_graph.WCGraph.pull_out_var_hecke"></a>

#### pull\_out\_var\_hecke

```python
@classmethod
def pull_out_var_hecke(cls, w: Permutation,
                       k: int) -> set[tuple[tuple[int, ...], Permutation]]
```

Hecke/Demazure analogue of pull_out_var(1, w) for fixed first-row size.

Returns all pairs ``(row, wpp)`` with:
- ``row`` a strictly decreasing tuple of simple reflections of size ``k``
- ``wpp`` a permutation such that
  ``Permutation.ref_product(*row) @ wpp.shiftup(1) == w``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.all_wc_graphs_slow"></a>

#### all\_wc\_graphs\_slow

```python
@classmethod
def all_wc_graphs_slow(cls,
                       perm: Permutation,
                       length: int = -1,
                       weight: tuple[int, ...] | None = None) -> set[WCGraph]
```

Generate all WC graphs with Demazure permutation ``perm`` and fixed row count.

The recursion mirrors ``all_rc_graphs`` but uses Hecke/Demazure pull-out on the
first row via ``pull_out_var_hecke``.

<a id="schubmult.combinatorics.wc_graph.WCGraph.grothendieck_polynomial_via_wc"></a>

#### grothendieck\_polynomial\_via\_wc

```python
@classmethod
def grothendieck_polynomial_via_wc(cls,
                                   perm: Permutation,
                                   x: Sequence[Expr],
                                   beta: Expr,
                                   length: int = -1)
```

Compute a Grothendieck candidate by summing WC graph monomials.

Uses weight monomials with ``beta`` exponent ``|word|-inv(perm)``.

<a id="schubmult.combinatorics.plactic"></a>

# schubmult.combinatorics.plactic

<a id="schubmult.combinatorics.plactic.Plactic"></a>

## Plactic Objects

```python
class Plactic(GridPrint, CrystalGraph)
```

<a id="schubmult.combinatorics.plactic.Plactic.up_jdt_slide"></a>

#### up\_jdt\_slide

```python
def up_jdt_slide(row, col)
```

Perform an upward jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new Plactic tableau.

<a id="schubmult.combinatorics.plactic.Plactic.down_jdt_slide"></a>

#### down\_jdt\_slide

```python
def down_jdt_slide(row, col)
```

Perform a jeu de taquin slide starting from the given (row, col)
position (0-indexed). Returns a new Plactic tableau.

<a id="schubmult.combinatorics.plactic.Plactic.shiftup"></a>

#### shiftup

```python
def shiftup(k)
```

Return a new Plactic with all entries increased by k.

<a id="schubmult.combinatorics.plactic.Plactic.all_ss_tableaux"></a>

#### all\_ss\_tableaux

```python
@classmethod
def all_ss_tableaux(cls, shape, max_entry, inner_shape=None)
```

Generate all semistandard tableaux of given shape (or skew shape) with entries <= max_entry.

**Arguments**:

  
- `shape` - Sequence of row lengths (outer shape)
- `max_entry` - Maximum entry value
- `inner_shape` - Optional sequence of left offsets per row (for skew shapes).
  If provided, positions [row][0:inner_shape[row]] are marked as 0.
  

**Returns**:

  Set of Plactic instances representing all valid semistandard tableaux

<a id="schubmult.combinatorics.plactic.Plactic.row_word"></a>

#### row\_word

```python
@property
def row_word()
```

Return the row-reading word as a flat tuple.

<a id="schubmult.combinatorics.plactic.Plactic.transpose"></a>

#### transpose

```python
def transpose()
```

Return the transpose of this Plactic tableau.

<a id="schubmult.combinatorics.plactic.Plactic.invert"></a>

#### invert

```python
def invert()
```

Return a Plactic whose entries are remapped so that standard
(increasing) insertion order applies. If reverse_semistandard is True
we negate entries (so larger original becomes smaller).

<a id="schubmult.combinatorics.plactic.Plactic.__mul__"></a>

#### \_\_mul\_\_

```python
def __mul__(other)
```

Plactic product: insert entries of `other` in row-reading order
(top-to-bottom, left-to-right) into a copy of self.

<a id="schubmult.combinatorics.plactic.Plactic.skew_shape"></a>

#### skew\_shape

```python
@property
def skew_shape()
```

Return the skew shape as a tuple of (row_length, left_offset) pairs.

<a id="schubmult.combinatorics.plactic.Plactic.rs_insert"></a>

#### rs\_insert

```python
def rs_insert(*letters)
```

Insert one or more letters in sequence (row-insertion) and return a new Plactic.

<a id="schubmult.combinatorics.plactic.Plactic.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

Crystal raising operator e_i on the Plactic tableau (delegates to RCGraph).

<a id="schubmult.combinatorics.plactic.Plactic.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i)
```

Crystal lowering operator f_i on the Plactic tableau (delegates to RCGraph).

<a id="schubmult.combinatorics.plactic.Plactic.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Return the crystal weight of this tableau (delegated to RCGraph).

<a id="schubmult.combinatorics.plactic.Plactic.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

Return the length/number of rows used for the crystal

<a id="schubmult.combinatorics.plactic.Plactic.yamanouchi"></a>

#### yamanouchi

```python
@classmethod
def yamanouchi(cls, shape)
```

Return the Yamanouchi (highest-weight) tableau of the given shape.

<a id="schubmult.combinatorics.plactic.Plactic.is_increasing"></a>

#### is\_increasing

```python
@property
def is_increasing()
```

Check if the tableau is strictly increasing in rows and columns.

<a id="schubmult.combinatorics.plactic.Plactic.reverse_rsk"></a>

#### reverse\_rsk

```python
def reverse_rsk(recording_tableau)
```

Inverse RSK (row-insertion) for the pair (P,Q) where `self` is P and
`recording_tableau` is the standard recording tableau Q of the same shape.

Returns the original word as a list of integers (in insertion order).

<a id="schubmult.combinatorics.plactic.Plactic.rsk_insert"></a>

#### rsk\_insert

```python
@classmethod
def rsk_insert(cls, *letters)
```

Perform ordinary RSK (row insertion) on the given sequence of letters,
starting from this Plactic as the initial P-tableau. Returns a pair
(P_tableau, Q_tableau) where both are Plactic instances and Q is the
standard recording tableau with entries 1..m (in insertion order).

Usage:
  P, Q = Plactic().rsk_insert(3,1,2,1)
  or
  P, Q = Plactic().rsk_insert([3,1,2,1])

<a id="schubmult.combinatorics.plactic.Plactic.reverse_rectify_to_outer"></a>

#### reverse\_rectify\_to\_outer

```python
def reverse_rectify_to_outer(outer_shape)
```

Deterministic reverse-rectification to a given outer shape `outer_shape`.

Given a (straight) tableau `self` of shape lambda, produce a skew tableau
(represented as a Plactic whose rows may contain 0's for inner cells)
of outer shape `outer_shape` whose rectification is `self`.

outer_shape: iterable of nonnegative ints giving the desired outer row
lengths (mu_0 >= mu_1 >= ...).

Algorithm (deterministic):

- Let mu be the set of cells (r,c) with 0 <= r < len(mu) and 0 <= c < mu[r].
- While the current set of occupied cells (from the working tableau) is
a strict subset of mu:

* choose an outer corner cell (r,c) in mu\current_cells (no cell of mu
to its right or below). Choose the maximal such (r,c) (deterministic).
* create a hole at (r,c) (extend rows/cols as needed, set that cell to 0),
then perform an upward jeu-de-taquin slide from (r,c) using
up_jdt_slide to move the hole inward.
* adopt the resulting tableau and continue.

- Return the resulting Plactic (with zeros marking inner/removed cells).

**Notes**:

  
  - Raises ValueError if outer_shape does not dominate the current shape
  (i.e. mu must contain the current occupied cells).
  - Raises RuntimeError if no suitable outer corner can be found or if an
  up_jdt_slide fails (this indicates the requested outer shape is not
  attainable by reverse-rectification).

<a id="schubmult.combinatorics.hecke_plactic"></a>

# schubmult.combinatorics.hecke\_plactic

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

<a id="schubmult.mult.single"></a>

# schubmult.mult.single

Ordinary (single) Schubert polynomial multiplication.

Implements the ``schubmult_py`` kernel: the product of a linear combination of
Schubert polynomials ``S_u`` (given as a dict ``{u: coeff}``) with a single
Schubert polynomial ``S_v``, returned as a coefficient dict ``{w: coeff}``.

The algorithm is the recursive "v-path" / transition method used throughout
``schubmult``: write ``theta = (~v).theta()`` for the dominant weakly-decreasing
vector bounding ``v``'s Lehmer code, let ``mu = uncode(theta)`` and
``vmu = v * mu``, and process the entries of ``theta`` one at a time, tracking
Bruhat-chain "v-paths" from ``vmu`` down to the identity (``compute_vpathdicts``)
alongside chains of elementary-symmetric moves on the ``u`` side
(``elem_sym_perms``). Accumulating consistent pairs of chains and reading off
the coefficient landing on ``vmu`` gives the product.

<a id="schubmult.mult.single.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum)
```

Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum``.

Uses the classical Monk rule: ``x_k * S_u = sum S_{u t_{ij}}`` over Bruhat
covers ``u t_{ij}`` with ``i <= k < j`` (added) minus those with ``j <= k < i``
(subtracted), via ``elem_sym_perms(u, 1, varnum)``.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `varnum` - 1-indexed variable index ``k``.
  

**Returns**:

- `dict` - The updated coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.mult_poly_py"></a>

#### mult\_poly\_py

```python
def mult_poly_py(coeff_dict, poly, var_x=_vars.var_x)
```

Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``; each single
variable leaf is dispatched to ``single_variable``, and any other leaf just
scales every coefficient.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `poly` - Symbolic polynomial expression in the variables of ``var_x``.
- `var_x` - Generating set identifying the ``x`` variables (default ``x``).
  

**Returns**:

- `dict` - The updated coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.schubmult_py"></a>

#### schubmult\_py

```python
def schubmult_py(perm_dict, v)
```

Multiply ``sum_u coeff_u S_u(x)`` by the (ordinary) Schubert polynomial ``S_v``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available and the
permutations fit within its ``MAXN``, falling back to the pure-Python
implementation otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}`` (integer coefficients).
- `v` - Permutation (or array-form list) indexing the Schubert polynomial to
  multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}`` for the product.

<a id="schubmult.mult.single.schubmult_py_down"></a>

#### schubmult\_py\_down

```python
def schubmult_py_down(perm_dict, v)
```

Divided-difference ("down") variant of ``_schubmult_py_python``.

Same v-path recursion but built from ``elem_sym_perms_op`` (Bruhat *descents*)
instead of ``elem_sym_perms``, used for the down/dual side of the transition
recursion rather than ordinary multiplication.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation to multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.single.schub_coprod_py"></a>

#### schub\_coprod\_py

```python
def schub_coprod_py(perm, indices)
```

Coproduct of ``S_perm`` restricted to the variable split named by ``indices``.

Computes the expansion of the (single) Schubert polynomial coproduct
``Delta_{indices}(S_perm) = sum (firstperm, secondperm) -> coeff`` by
multiplying the Grassmannian permutation for ``indices`` against ``perm``
via ``schubmult_py`` and splitting each resulting permutation's window into
its first ``N`` and remaining ``len(perm) - N`` values.

**Arguments**:

- `perm` - Permutation (or array-form list) to take the coproduct of.
- `indices` - Iterable of 1-indexed positions selecting the variable split.
  

**Returns**:

- `dict` - Mapping ``{(firstperm, secondperm): coeff}``.

<a id="schubmult.mult.separated_descents"></a>

# schubmult.mult.separated\_descents

Separated-descents product of (double) Grothendieck polynomials.

Implements the *pipe puzzle* formula of Fan--Guo--Xiong,
"Bumpless pipe dreams meet puzzles" (arXiv:2309.00467), Theorem 2.5.

For permutations ``u`` and ``v`` with *separated descents* at a position ``k``
(i.e. every descent of ``u`` is ``<= k`` and every descent of ``v`` is ``>= k``)
the double Grothendieck polynomials satisfy

    G_u(x, y) * G_v(x, t) = sum_w c_{u,v}^w(t, y) * G_w(x, t),

and the structure constants ``c_{u,v}^w(t, y)`` are given by a positive
(in the sense of Theorem 2.5) sum over pipe puzzles.

The beta convention matches the rest of ``schubmult``: the multiplicative
formal group law is ``x (+) y = x + y + beta*x*y`` with formal subtraction

    x (-) y = (x - y) / (1 + beta*y).

Setting ``beta = 0`` recovers the double Schubert (cohomology) structure
constants of Theorem 4.x (the "Schubert pipe puzzle" specialization).

<a id="schubmult.mult.separated_descents.separated_descents_coeffs"></a>

#### separated\_descents\_coeffs

```python
def separated_descents_coeffs(u, v, var1, var2, beta=None, grid_size=None)
```

Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
``var2`` are the ``y`` variables (secondary variables of ``u``).

Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

In K-theory the product may involve ``G_w`` with ``w`` in a larger
symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
``grid_size`` is not supplied it is chosen large enough to capture every
such ``w``: the maximum value of any appearing ``w`` is bounded by
``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.

<a id="schubmult.mult.separated_descents.separated_descents_coeffs_plus"></a>

#### separated\_descents\_coeffs\_plus

```python
def separated_descents_coeffs_plus(u,
                                   v,
                                   var1,
                                   var2,
                                   beta=None,
                                   grid_size=None,
                                   mangle_genset=False)
```

Coefficients ``c_{u,v}^w(var1, var2)`` for a single pair ``u``, ``v``.

``var1`` are the ``t`` variables (secondary variables of ``v`` and ``w``),
``var2`` are the ``y`` variables (secondary variables of ``u``).

Returns a dict ``{w: c_{u,v}^w}`` with symbolic coefficients.

In K-theory the product may involve ``G_w`` with ``w`` in a larger
symmetric group ``S_{n'}`` (Remark following Theorem 2.5).  When
``grid_size`` is not supplied it is chosen large enough to capture every
such ``w``: the maximum value of any appearing ``w`` is bounded by
``(n - 1) + deg(G_u) + deg(G_v)`` where ``n = max(len(u), len(v))``.

<a id="schubmult.mult.separated_descents.grothmult_double"></a>

#### grothmult\_double

```python
def grothmult_double(perm_dict, v, var1, var2, beta=None)
```

Separated-descents product of double Grothendieck polynomials.

Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

    coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.

<a id="schubmult.mult.separated_descents.grothmult_double_plus"></a>

#### grothmult\_double\_plus

```python
def grothmult_double_plus(perm_dict,
                          v,
                          var1,
                          var2,
                          beta=None,
                          mangle_genset=False)
```

Separated-descents product of double Grothendieck polynomials.

Given ``perm_dict = {u: coeff_u}`` and a permutation ``v`` such that every
``u`` has separated descents with ``v``, returns ``{w: coeff_w}`` where

    coeff_w = sum_u c_{u,v}^w(var1, var2) * coeff_u,

with ``c_{u,v}^w`` the pipe-puzzle structure constants of Theorem 2.5.

``var1`` are the ``t`` variables, ``var2`` the ``y`` variables.

<a id="schubmult.mult.groth"></a>

# schubmult.mult.groth

beta-Grothendieck Chevalley formula: multiplication of a (single) Grothendieck
polynomial by a bare x_k variable, i.e. x_k * G_w^(beta).

Non-equivariant (y=0) for now; ``GrothendieckRing`` has no coefficient/y genset yet.

Derived from M. Willems, "A Chevalley formula in equivariant K-theory"
(arXiv:math/0603220), Theorem 5 (the ordinary, non-equivariant specialization of
his equivariant Chevalley formula, Theorem 4). Willems indexes K-theory classes
O_w by the *dimension* of the Schubert variety, dual to the *codimension*
indexing used by Schubert/Grothendieck polynomials S_w/G_w; the w0-conjugation
below (``hat_w = w0*w`` going in, ``w0*v`` coming out) translates between the two
conventions. The beta-grading (beta^(d-1) per length difference d = l(v)-l(w))
matches this codebase's beta-deformed Grothendieck polynomial normalization
(beta=0 recovers the classical double Schubert Monk formula). Calibrated against
grothendieck_poly()/to_groth() (see session notes).

<a id="schubmult.mult.groth.chevalley_x_k"></a>

#### chevalley\_x\_k

```python
def chevalley_x_k(w, k, beta, n=None)
```

Coefficients of ``x_k * G_w^(beta)`` in the Grothendieck basis, as a dict
``{v: coeff}`` (``w`` itself never appears: the self-term cancels identically).

<a id="schubmult.mult.groth.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by the single variable ``x_varnum``
(Grothendieck Chevalley formula), via ``chevalley_x_k``.

<a id="schubmult.mult.groth.mult_poly_groth"></a>

#### mult\_poly\_groth

```python
def mult_poly_groth(coeff_dict, poly, var_x, beta)
```

Multiply ``sum_u coeff_u G_u^(beta)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable_groth``.

<a id="schubmult.mult.positivity"></a>

# schubmult.mult.positivity

Manifestly positive representations of double Schubert structure constants.

The structure constants ``c^w_{u,v}(y, z)`` of double Schubert polynomial
multiplication (``schubmult_double``) are known to be polynomials in the
differences ``y_i - z_j`` with nonnegative integer coefficients (Graham's
positivity theorem). This module computes that manifestly positive form:

- ``posify``: the main recursive engine. Reduces ``(u, v, w)`` via known
  combinatorial identities (pattern-avoidance checks, dominance/one-dominance,
  descent and coefficient reductions in ``schubmult.utils.schub_lib``) down to
  cases with closed positive formulas (``dualcoeff``, ``forwardcoeff``, or a
  single elementary symmetric polynomial), falling back to the integer-LP
  solver ``compute_positive_rep`` when no reduction applies.
- ``compute_positive_rep``: expresses an arbitrary such polynomial as a
  nonnegative-integer combination of product-of-differences monomials, found
  via an integer program (PuLP) over a spanning set of candidate monomials.
- ``dualcoeff``/``forwardcoeff``/``dualpieri``: closed-form positive rules for
  special cases (``u`` dominates ``w``, the ``will_formula_work`` forward Monk
  case, and the dual Pieri expansion respectively).

<a id="schubmult.mult.positivity.compute_positive_rep"></a>

#### compute\_positive\_rep

```python
def compute_positive_rep(val, var2=None, var3=None, msg=False)
```

Express ``val`` as a nonnegative-integer combination of product-of-differences monomials.

``val`` must be a polynomial in ``var2``/``var3`` known (by positivity of
double Schubert structure constants) to admit an expansion
``sum_b n_b * prod (var2_i - var3_j)`` with ``n_b >= 0`` integers. Builds a
candidate spanning set of such product monomials from ``val``'s own
monomials, then solves an integer program (via PuLP) for nonnegative
integer coefficients ``n_b`` matching ``val`` exactly.

**Arguments**:

- `val` - Symbolic polynomial expression in ``var2``/``var3``.
- `var2` - First secondary alphabet (``y``).
- `var3` - Second secondary alphabet (``z``).
- `msg` - Passed through to the LP solver as its ``msg`` (verbosity) option.
  

**Returns**:

  A symbolic expression equal to ``val``, written as a sum of
  nonnegative-integer multiples of product-of-differences monomials.
  

**Raises**:

- `Exception` - If the reconstructed expression does not equal ``val``
  (i.e. no valid nonnegative integer solution reproduces it exactly).

<a id="schubmult.mult.positivity.posify"></a>

#### posify

```python
@cached(
    cache={},
    key=lambda val, u2, v2, w2, var2=None, var3=
    None, msg=False, sign_only=False, optimize=True: hashkey(
        val, u2, v2, w2, var2, var3, msg, sign_only, optimize),
)
def posify(val,
           u2,
           v2,
           w2,
           var2=None,
           var3=None,
           msg=False,
           sign_only=False,
           optimize=True,
           n=_vars.n)
```

Manifestly positive representation of the structure constant ``c^{w2}_{u2,v2}(var2, var3)``.

``val`` is the (already computed, possibly not manifestly positive) value of
the coefficient of ``S_{w2}`` in ``S_{u2}(x, var2) * S_{v2}(x, var3)``.
Recursively reduces ``(u2, v2, w2)`` via pattern-avoidance-guarded identities
(``try_reduce_u``/``try_reduce_v``, ``reduce_descents``, ``reduce_coeff``,
``is_split_two``) toward cases handled by closed positive formulas
(a single elementary symmetric polynomial when ``v`` has one nonzero code
entry, ``dualcoeff`` when ``will_formula_work(v, u)`` or ``u`` dominates
``w``, ``forwardcoeff`` when ``will_formula_work(u, v)``, or the
length-one-difference case built from ``pull_out_var``/``schubpoly``
directly). Falls back to ``compute_positive_rep`` (an integer-LP search)
when no reduction or closed formula applies and ``optimize`` is true.

Results are cached by ``(val, u2, v2, w2, var2, var3, msg, sign_only, optimize)``.

**Arguments**:

- `val` - The structure constant to re-express positively.
- `u2` - First factor's permutation.
- `v2` - Second factor's permutation.
- `w2` - Target permutation (coefficient of ``S_{w2}``).
- `var2` - First secondary alphabet.
- `var3` - Second secondary alphabet.
- `msg` - Verbosity flag passed down to ``compute_positive_rep``'s LP solver.
- `sign_only` - If ``True``, only determine and return the sign of ``val``
  (``-1``, ``0``, or ``1``) rather than a full positive expression.
- `optimize` - If ``False``, return ``val`` unchanged when no closed-form
  reduction applies (skip the LP fallback); if ``None`` and that
  case is reached, raise.
- `n` - Size of the ambient alphabet used when no other bound is available.
  

**Returns**:

  A manifestly positive expression equal to ``val`` (or, if
  ``sign_only``, one of ``-1``, ``0``, ``1``).

<a id="schubmult.mult.positivity.shiftsub"></a>

#### shiftsub

```python
def shiftsub(pol, var2=None)
```

Shift every ``var2[i]`` in ``pol`` up to ``var2[i + 1]`` (for ``i`` in ``0..98``).

<a id="schubmult.mult.positivity.posify_generic_partial"></a>

#### posify\_generic\_partial

```python
def posify_generic_partial(val, u2, v2, w2)
```

``posify`` specialized to the fixed generic alphabets ``_vars.var_g1``/``_vars.var_g2``.

Asserts (raises on mismatch) that the recomputed positive expression equals
the input ``val``, as a consistency check.

<a id="schubmult.mult.positivity.schubmult_generic_partial_posify"></a>

#### schubmult\_generic\_partial\_posify

```python
@cache
def schubmult_generic_partial_posify(u2, v2)
```

Manifestly positive expansion of ``S_{u2}(x, var_g1) * S_{v2}(x, var_g2)``.

Returns ``{w2: coeff}`` where each ``coeff`` is the positive representation
(via ``posify_generic_partial``) of the corresponding
``schubmult_double_pair_generic_alt`` coefficient.

<a id="schubmult.mult.positivity.forwardcoeff"></a>

#### forwardcoeff

```python
def forwardcoeff(u, v, perm, var2=None, var3=None)
```

Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "forward" case
(used when ``will_formula_work(u, v)`` holds in ``posify``).

Writes ``muv = uncode(v.theta())`` and reduces to a lookup in
``schubmult_double_pair(u, muv, var2, var3)`` when the length condition
``(perm * (~v * muv)).inv == (~v * muv).inv + perm.inv`` holds; returns 0 otherwise.

<a id="schubmult.mult.positivity.dualcoeff"></a>

#### dualcoeff

```python
def dualcoeff(u, v, perm, var2=None, var3=None)
```

Closed-form structure constant ``c^{perm}_{u,v}(var2, var3)`` for the "dual" case
(used in ``posify`` when ``will_formula_work(v, u)`` holds or ``u`` dominates ``perm``).

When ``u`` is the identity, reduces directly to a single Schubert
polynomial ``schubpoly(v * (~perm), var2, var3)``. Otherwise expands via
``dualpieri`` (directly if ``u`` dominates ``perm``, or after rewriting to
``u``'s dominant permutation ``uncode(u.theta())`` otherwise), summing
products of ``(var2_{i+1} - var3_j)`` factors against a final
``schubpoly`` term.

<a id="schubmult.mult.positivity.dualpieri"></a>

#### dualpieri

```python
def dualpieri(mu, v, w)
```

Dual Pieri expansion used by ``dualcoeff``: enumerate the data witnessing
``S_mu * S_v -> S_w`` when ``mu`` is dominant.

Compares ``mu``'s inverse code against ``w``'s inverse code layer by layer,
peeling one "cycle" of variables per layer via ``divdiffable``/``pull_out_var``,
and returns the list of ``[vlist, vp]`` pairs consumed by ``dualcoeff`` to
build the final positive expression (empty list if ``w`` is not reachable
from ``mu``, ``v`` this way).

<a id="schubmult.mult.quantum"></a>

# schubmult.mult.quantum

Quantum (single) Schubert polynomial multiplication.

Implements ``schubmult_q``/``schubmult_q_fast``: the product of a linear
combination of quantum Schubert polynomials ``S_u`` with a single ``S_v``,
returned as a coefficient dict ``{w: coeff}`` polynomial in the quantum
parameters ``q_1, q_2, ...``. Uses the same ``theta``/v-path recursion as
``schubmult.mult.single``, with the elementary-symmetric step generalized to
the quantum Pieri-type moves of ``elem_sym_perms_q`` (which may pick up a
factor of ``q`` when a Bruhat move is replaced by its quantum analogue).
``schubmult_q`` uses ``strict_theta`` (no repeated layer merging);
``schubmult_q_fast``/``_schubmult_q_fast_python`` uses ``medium_theta`` and
merges adjacent equal-length layers via ``double_elem_sym_q`` for speed.

<a id="schubmult.mult.quantum.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var_q=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by the single variable ``x_varnum`` (quantum Monk rule).

Same structure as ``schubmult.mult.single.single_variable``, using
``elem_sym_perms_q`` so that some Bruhat moves carry a factor from ``var_q``.

<a id="schubmult.mult.quantum.mult_poly_q"></a>

#### mult\_poly\_q

```python
def mult_poly_q(coeff_dict, poly, var_x=_vars.var_x, var_q=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_py``.

<a id="schubmult.mult.quantum.schubmult_q_fast"></a>

#### schubmult\_q\_fast

```python
def schubmult_q_fast(perm_dict, v, q_var=_vars.q_var)
```

Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available, falling
back to ``_schubmult_q_fast_python`` (the ``medium_theta``-based recursion
with merged equal-length layers) otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the quantum Schubert
  polynomial to multiply by.
- `q_var` - Generating set for the quantum parameters.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``, polynomial in ``q_var``.

<a id="schubmult.mult.quantum.schubmult_q"></a>

#### schubmult\_q

```python
def schubmult_q(perm_dict, v)
```

Multiply ``sum_u coeff_u S_u(x)`` by the quantum Schubert polynomial ``S_v``.

Reference (non-"fast") implementation: uses ``strict_theta`` and processes
every layer individually (no merging of equal-length adjacent layers), so it
is simpler but slower than ``schubmult_q_fast``. Results agree with
``schubmult_q_fast`` for all inputs.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the quantum Schubert
  polynomial to multiply by.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``, polynomial in the
  default quantum parameters ``q``.

<a id="schubmult.mult.groth_double"></a>

# schubmult.mult.groth\_double

K-theoretic Monk formula for double Grothendieck polynomials.

Implements Lenart--Postnikov, *Affine Weyl groups in K-theory and representation
theory* (arXiv:math/0309207), Corollary 8.2 (the :math:`K_T`-Monk formula), in
type :math:`A` and transported to the ``schubmult`` conventions.

Reconciling the conventions
---------------------------
``schubmult`` has no exponentials: it uses the multiplicative formal group law
``a (+) b = a + b + beta*a*b`` with formal inverse ``(-)b = -b/(1 + beta*b)``.
The double Grothendieck polynomial of a simple reflection is

.. math::

    \mathfrak{G}_{s_k}(x, y)
        = \frac{1}{\beta}\Bigl(\prod_{i=1}^{k}(1 + \beta x_i)(1 + \beta y_i) - 1\Bigr),

which at ``beta = -1`` is precisely the class ``1 - x^{w_0(omega_k)} e^{-omega_k}``
of Lemma 8.1(a).  The dictionary between the paper and this package is

* torus characters: ``x^{eps_i}  <->  1 + beta*y_i`` (so that ``(x^{eps_a - eps_b} - 1)/beta``
  is the root ``y_a (-) y_b``, matching ``DoubleGrothendieckRing.exp_root``);
* basis elements: ``[O_{X_{w_0 w}}]  <->  (-beta)^{l(w)} G_w``, since the paper
  indexes structure sheaves by the *dimension* of the Schubert variety whereas
  ``G_w`` has lowest term ``S_w`` of degree ``l(w)`` (codimension indexing).

Under ``u -> w_0 u`` the saturated *decreasing* chains of Corollary 8.2 become
saturated *increasing* chains, the reflections ``t_{ij}`` are unchanged, and the
character prefactor ``x^{nu(J)} = x^{w_0(omega_k) - u(omega_k)}`` (constant in
``J`` because ``omega_k`` is minuscule in type ``A``) becomes
``prod_{i<=k} (1 + beta*y_i)/(1 + beta*y_{u(i)})``.  Rescaling by
``(-beta)^{l(.)}`` turns the signs ``(-1)^{|J|}`` into powers ``beta^{|J|}`` and
yields

.. math::

    \mathfrak{G}_u(x, y)\,\mathfrak{G}_{s_k}(x, z)
        = \frac{1}{\beta}\Bigl(\Theta_u \sum_J \beta^{|J|}\,
          \mathfrak{G}_{u\,r_J}(x, y) - \mathfrak{G}_u(x, y)\Bigr),
    \qquad
    \Theta_u = \prod_{i=1}^{k}\frac{1 + \beta z_i}{1 + \beta y_{u(i)}},

the sum being over the subsets ``J`` of a reduced ``(-omega_k)``-chain of
reflections whose reflections build a saturated increasing Bruhat chain from
``u`` (the empty subset included).  The two secondary alphabets are handled by
``G_{s_k}(x, z) = C G_{s_k}(x, y) + (C - 1)/beta`` with
``C = prod_{i<=k}(1 + beta*z_i)/(1 + beta*y_i)``, which is exactly what turns
the ``y_i`` of ``x^{nu}`` into the ``z_i`` of ``Theta_u``.

By Corollary 15.4 of the same paper a reduced ``(-omega_k)``-chain of
reflections in type ``A_{n-1}`` is

    ``t_{1,n}, t_{1,n-1}, ..., t_{1,k+1}, t_{2,n}, ..., t_{2,k+1}, ..., t_{k,k+1}``.

<a id="schubmult.mult.groth_double.monk_chain"></a>

#### monk\_chain

```python
def monk_chain(k)
```

Reduced ``(-omega_k)``-chain of reflections in ``A_{n-1}``, as ``(i, j)``, ``i <= k < j``.

``omega_k = eps_1 + ... + eps_k``, so this is ``epsilon_chain`` on ``{1, ..., k}``:
every root ``alpha_{ij}`` with ``i, j <= k`` pairs to zero and drops out, and the
survivors all sit at level ``1``.  The order is ``i`` decreasing, then ``j`` decreasing.

<a id="schubmult.mult.groth_double.epsilon_chain"></a>

#### epsilon\_chain

```python
def epsilon_chain(positions, inverse=False, ambient_rank=None)
```

Reduced ``(-eps_A)``-chain of reflections in ``A_{n-1}``, ``A = positions``.

``positions`` is a single index or an iterable of them (repeats allowed), and
``eps_A = sum_{i in A} eps_i``.  With ``inverse=True`` the chain is for ``+eps_A``
instead, which is the weight of the inverse class ``prod_{i in A}(1 + beta*x_i)^{-1}``.

``(-omega_k)``-chains only see the roots ``alpha_{ij}`` with ``i <= k < j``, which is
why ``monk_chain`` multiplies by the whole product ``prod_{i<=k}(1 + beta*x_i)``.
Selecting an arbitrary set of variables needs ``eps_A`` instead, whose chain also
involves the roots ``alpha_{ik}`` with ``i < k``.

Built by Prop. 6.7: the reflections ``s_{alpha, m}`` separating the fundamental
alcove from ``A_{eps_A}``, ordered by the lexicographic key
``(lambda, alpha^vee)^{-1} (-m, (omega_1, alpha^vee), ..., (omega_{n-1}, alpha^vee))``.
Entries are ``(a, b, m)`` for the positive root ``alpha_{ab} = eps_a - eps_b``, ``a < b``;
``m > 0`` means ``b(r) = -alpha`` is negative and the step carries a sign in Thm 6.1.

Concatenating the individual ``(-eps_i)``-chains would also be legal (Prop. 12.2) but
only after translating the blocks, which shifts their levels; going through Prop. 6.7
avoids that and is reduced.  Note a single ``k`` gives ``(i, k, 0)`` for ``i < k`` and
``(k, j, 1)`` for ``j > k``, while for ``A = {1, ..., k}`` every ``alpha_{ij}`` with
``i, j <= k`` pairs to zero and drops out, leaving exactly ``monk_chain(k)``.
Flipping to ``inverse=True`` exchanges those two families.

<a id="schubmult.mult.groth_double.one_plus_beta_x_groth"></a>

#### one\_plus\_beta\_x\_groth

```python
def one_plus_beta_x_groth(coeff_dict,
                          positions,
                          var2=None,
                          beta=None,
                          inverse=False)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by ``prod_{i in positions} (1 + beta*x_i)``.

The one-pass Pieri rule of Theorem 6.1 at ``lambda = -eps_A``; see
``_one_plus_beta_x_terms`` for the coefficient.  ``inverse=True`` gives the inverse
operator ``prod_{i in positions} (1 + beta*x_i)^{-1}``, i.e. ``lambda = +eps_A``.

<a id="schubmult.mult.groth_double.single_variable_groth"></a>

#### single\_variable\_groth

```python
def single_variable_groth(coeff_dict, varnum, var2=None, beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by the single variable ``x_varnum``.

Returns ``{w: coeff_w}``.  This is ``_one_plus_beta_x_terms`` with the
identity subtracted off and ``beta`` divided out; the diagonal coefficient
collapses to ``(-) var2[u(varnum)] = -y/(1 + beta*y)``, the formal inverse of
``var2[u(varnum)]``, which is the localization of ``x_varnum`` at ``u``.

<a id="schubmult.mult.groth_double.elem_sym_perms_groth"></a>

#### elem\_sym\_perms\_groth

```python
def elem_sym_perms_groth(u, k)
```

K-theoretic analogue of ``elem_sym_perms``: ``{w: {d: multiplicity}}``.

Same recursion as ``elem_sym_perms(u, p, k)`` -- a step is any Bruhat cover
``w -> w t_{ij}`` with ``i <= k < j``, and ``j`` is required to weakly decrease along
the chain -- but with the ``p`` cut-off dropped, so chains of every length are
produced and the degree cut-off is left to the coefficient.

This is deliberately *not* a subset-of-a-fixed-chain enumeration.  A
``lambda``-chain imposes a total order on the transpositions, which loses covers such
as ``id < [1,3,2] < [2,3,1]`` (that needs ``t_{23}`` before ``t_{13}``); covers come
from arbitrary upward transpositions, and only ``j`` is constrained.

``d = l(w) - l(u)`` is the chain length.  A position ``i <= k`` may be stepped on more
than once, and distinct chains can land on the same ``w`` at the same ``d``, which is
the source of the K-theoretic multiplicities.

<a id="schubmult.mult.groth_double.grothmult_double_block"></a>

#### grothmult\_double\_block

```python
def grothmult_double_block(coeff_dict,
                           positions,
                           zvar=None,
                           var2=None,
                           beta=None,
                           fgl=True)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by a linear block over ``positions``:

    fgl=True  ->  prod_{i in A} (x_i (+) zvar),   x (+) z = x*(1 + beta*z) + z
    fgl=False ->  prod_{i in A} (x_i - zvar)

``positions`` is an arbitrary index set (repeats allowed), matching the ``index_list``
that ``pull_out_var`` produces, so this is the ``G``-basis analogue of the top-degree
mixed-variable block driving ``schubmult_double_alt`` / ``DoubleSchubertRing.elem_mul``.

``fgl=False`` is the plain ``beta = 0`` block ``(x_1 - z)(x_2 - z)...`` -- still a
perfectly good operator on the ``G`` basis, and the two are interchangeable via
``x (+) z = (1 + beta*z) * (x - (-)z)``, so either can be recovered from the other by
rescaling ``zvar``.

Computed by folding ``single_variable_groth`` one position at a time, using
``(a x_i + b) F = a (x_i F) + b F``.  That keeps every intermediate coefficient
polynomial in ``beta``; the one-pass alternative would expand
``prod_i ((1 + beta*x_i)(1 + beta*z) - 1) / beta**|A|`` by inclusion-exclusion over the
subsets of ``A`` (each term a ``one_plus_beta_x_groth`` call) and only cancel the
``beta^{-|A|}`` at the very end.

<a id="schubmult.mult.groth_double.groth_elem_sym_poly"></a>

#### groth\_elem\_sym\_poly

```python
def groth_elem_sym_poly(p, k, zvar, var_x, beta)
```

``E_p^beta(x_1..x_k; z) = e_p(x_1 (+) z, ..., x_k (+) z)``, ``x (+) z = x(1 + beta*z) + z``.

The double Grothendieck elementary symmetric: ``p == k`` gives
``(x_1(1 + beta*z) + z) ... (x_k(1 + beta*z) + z)`` and ``beta == 0`` gives the
factorial elementary symmetric ``elem_sym_poly(p, k, x, [-z])``.

<a id="schubmult.mult.groth_double.grothmult_double_pieri"></a>

#### grothmult\_double\_pieri

```python
def grothmult_double_pieri(coeff_dict,
                           p,
                           k,
                           zvar=None,
                           var_x=None,
                           var2=None,
                           beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by ``groth_elem_sym_poly(p, k, zvar, var_x, beta)``.

Exact, by folding the ``e_p`` DP over coefficient dicts one shifted variable
``x_i (+) zvar`` at a time (``single_variable_groth`` per step) -- no symbolic
expansion.  A closed-form Pieri rule in the style of ``dom_groth`` -- paths from
``elem_sym_perms_groth`` plus an ``elem_sym_poly`` in the localizations -- is *not*
implemented: grading the paths by ``beta^{d - m}`` with ``m`` the number of moved
positions and taking ``elem_sym_poly`` over the untouched ones is wrong already at
``p == k``.  The non-equivariant rule ``groth_pieri_mul`` grades instead by
``beta^{d - (number of marked steps)}`` with the multiplicity counting admissible
markings of the chain (``elem_sym_chains_groth``), so the equivariant coefficient
presumably needs that marking data rather than the moved/untouched split.

``var_x`` is unused (kept for signature compatibility).

<a id="schubmult.mult.groth_double.grothmult_double_top"></a>

#### grothmult\_double\_top

```python
def grothmult_double_top(coeff_dict, k, zvar=None, var2=None, beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var2)`` by the top linear block ``prod_{i=1}^{k}(x_i + zvar)``.

Closed positive Molev--Sagan Pieri rule (conjectural; verified exhaustively on
``S_4`` and sampled through ``S_6``, ``k <= 5``, against
``grothmult_double_block(..., zvar=-zvar, fgl=False)``):

.. math::

    \prod_{i=1}^{k}(x_i + z)\,\mathfrak{G}_u(x; y)
        = \sum_{w} \beta^{\,d - k + |Q|} \Bigl(\prod_{i=1}^{k} f_i\Bigr)\,
          \mathfrak{G}_w(x; y),
    \qquad d = \ell(w) - \ell(u),

where the factor ``f_i`` depends on the fate of the window value ``u(i)``:

* ``u(i) = w(i)`` (the set ``Q``):  ``(z(1 + beta*y_{u(i)}) - y_{u(i)}) / (1 + beta*y_{u(i)})``,
  i.e. ``z (+) (-)y_{u(i)}``, the K-theoretic analogue of ``z - y_{u(i)}``;
* ``u(i)`` stays in the window but moves left:  ``1 - beta*zvar``;
* ``u(i)`` exits the window or moves right within it:  ``1/(1 + beta*y_{u(i)})``.

The sum runs over the marked-chain K-Pieri support (``_top_block_support``).
At ``beta = 0`` this collapses to the ``p = k`` Pieri formula for double
Schubert polynomials [Samuel, Theorem 7.1]:
``S_u(x;y) prod(x_i - z) = sum_{u ->_k w} prod_{i in Q}(y_{u(i)} - z) S_w(x;y)``
with ``z -> -z``.

<a id="schubmult.mult.groth_double.mult_poly_groth_double"></a>

#### mult\_poly\_groth\_double

```python
def mult_poly_groth_double(coeff_dict,
                           poly,
                           var_x=None,
                           var_y=None,
                           beta=None)
```

Multiply ``sum_u coeff_u G_u(x, var_y)`` by an arbitrary polynomial in ``var_x``.

Mirrors ``mult_poly_double``; the leaves of the ``Add``/``Mul``/``Pow`` recursion
are handled by ``single_variable_groth``.

<a id="schubmult.mult.groth_double.dgroth_to_dschub"></a>

#### dgroth\_to\_dschub

```python
def dgroth_to_dschub(v, var3, beta=None)
```

Expand ``G_v(x, var3)`` in double Schubert polynomials: ``{v': coeff}``.

``sum_{v'} coeff_{v'} S_{v'}(x, var3) = G_v(x, var3)`` with coefficients in
``var3`` and ``beta``.  Exact but slow; delegates to ``grothendieck_poly``
with ``keep_as_schub=True``.

<a id="schubmult.mult.groth_double.groth_elem_sym_func"></a>

#### groth\_elem\_sym\_func

```python
def groth_elem_sym_func(k, i, u1, u2, v1, v2, vdiff, varl1, varl2, beta)
```

Expression form of ``_groth_elem_sym_frac``; see there for the rule.

<a id="schubmult.mult.groth_double.grothmult_double"></a>

#### grothmult\_double

```python
def grothmult_double(perm_dict, v, var2=None, var3=None, beta=None)
```

Multiply double Grothendieck polynomials, mirroring ``schubmult_double``.

Computes the expansion of ``sum_u coeff_u G_u(x, var2) * G_v(x, var3)`` in
the basis ``{G_w(x, var2)}`` and returns it as ``{w: coeff_w}``.

``v = s_k`` uses the verified chain formula of Corollary 8.2, and
``max_descent == 1`` folds that column by column.  General ``v`` goes through
``dgroth_to_dschub`` (exact, slow) and the conjectural vpath kernel
``_groth_schub_vpath_mul``, one run per double Schubert ``S_{v'}`` in the
expansion of ``G_v``.

The chain rank is inferred from the current permutation and selected positions.

<a id="schubmult.mult.double"></a>

# schubmult.mult.double

Double Schubert polynomial multiplication.

Implements the ``schubmult_double`` kernel: the product of a linear
combination of double Schubert polynomials ``S_u(x, var2)`` with a single
``S_v(x, var3)``, returned as a coefficient dict ``{w: coeff}`` of polynomials
in ``var2``/``var3``. Uses the same ``theta``/``vmu``/v-path recursion as
``schubmult.mult.single`` (see that module), replacing the plain elementary
symmetric contribution with ``elem_sym_func``, which carries the secondary
variables ``var2``/``var3``.

Also provides the "alt"/"from_elems" variants (building the product one
descent-pulled variable at a time via ``pull_out_var``, generic to any choice
of elementary-symmetric-like function), ``nilhecke_mult`` (nilHecke ring
multiplication), and ``schub_coprod_double`` (the double Schubert coproduct).

<a id="schubmult.mult.double.count_sorted"></a>

#### count\_sorted

```python
def count_sorted(mn, tp)
```

Count occurrences of ``tp`` in the sorted sequence ``mn`` via binary search.

<a id="schubmult.mult.double.single_variable"></a>

#### single\_variable

```python
def single_variable(coeff_dict, varnum, var2=None)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the single variable ``x_varnum``.

Equivariant Monk rule: the diagonal term contributes ``var2[u(varnum)]``
(localization of ``x_varnum`` at ``u``) and the off-diagonal terms are the
same Bruhat-cover moves as the ordinary (non-equivariant) ``single_variable``
in ``schubmult.mult.single``.

**Arguments**:

- `coeff_dict` - Mapping ``{Permutation: coeff}``.
- `varnum` - 1-indexed variable index ``k``.
- `var2` - Secondary (``y``) generating set.
  

**Returns**:

- `dict` - The updated coefficient dict.

<a id="schubmult.mult.double.single_variable_down"></a>

#### single\_variable\_down

```python
def single_variable_down(coeff_dict, varnum, var2=None)
```

Down (descent) variant of ``single_variable``, using ``elem_sym_perms_op``.

<a id="schubmult.mult.double.mult_poly_double"></a>

#### mult\_poly\_double

```python
def mult_poly_double(coeff_dict, poly, var_x=None, var_y=None)
```

Multiply ``sum_u coeff_u S_u(x, var_y)`` by an arbitrary polynomial ``poly`` in ``var_x``.

Recurses over the ``Add``/``Mul``/``Pow`` structure of ``poly``, dispatching
single-variable leaves to ``single_variable``; mirrors ``mult_poly_py`` with
the extra secondary alphabet ``var_y``.

<a id="schubmult.mult.double.mult_poly_double_alt"></a>

#### mult\_poly\_double\_alt

```python
def mult_poly_double_alt(coeff_dict, poly, var_x=None, var_y=None)
```

Variant of ``mult_poly_double`` that folds each factor via ``schubmult_double_dict``
instead of ``single_variable``, so ``poly`` is only ever expanded one variable/factor
at a time in the ``S_v`` basis rather than left as a raw scalar multiplier.

<a id="schubmult.mult.double.mult_poly_down"></a>

#### mult\_poly\_down

```python
def mult_poly_down(coeff_dict, poly)
```

Down (descent) variant of ``mult_poly_double``, using ``single_variable_down``
and the fixed default alphabet ``_vars.var1``.

<a id="schubmult.mult.double.nilhecke_mult"></a>

#### nilhecke\_mult

```python
def nilhecke_mult(coeff_dict1, coeff_dict2)
```

NilHecke ring product of ``coeff_dict1`` (polynomial coefficients) and
``coeff_dict2`` (permutation coefficients acting as divided-difference operators).

For each ``w`` in ``coeff_dict2`` its coefficient polynomial is pushed through
``mult_poly_down`` against ``coeff_dict1``, and each resulting permutation ``v``
is right-multiplied by ``w`` whenever that multiplication is length-additive.

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}``.

<a id="schubmult.mult.double.schubmult_double_pair"></a>

#### schubmult\_double\_pair

```python
@cache
def schubmult_double_pair(perm1, perm2, var2=None, var3=None)
```

``schubmult_double`` specialized to a single ``perm1`` with coefficient 1, cached.

<a id="schubmult.mult.double.schubmult_double_pair_generic"></a>

#### schubmult\_double\_pair\_generic

```python
@cache
def schubmult_double_pair_generic(perm1, perm2)
```

``schubmult_double_pair`` with the fixed generic secondary alphabets ``_vars.var_g1``/``_vars.var_g2``.

<a id="schubmult.mult.double.schubmult_double_pair_generic_alt"></a>

#### schubmult\_double\_pair\_generic\_alt

```python
@cache
def schubmult_double_pair_generic_alt(perm1, perm2)
```

Like ``schubmult_double_pair_generic`` but computed via ``schubmult_double_alt_from_elems``
with the factorial elementary symmetric function, then expanded/simplified.

<a id="schubmult.mult.double.schubmult_double_dict"></a>

#### schubmult\_double\_dict

```python
def schubmult_double_dict(perm_dict1, perm_dict2, var2=None, var3=None)
```

Multiply two coefficient dicts of double Schubert polynomials together.

Computes ``(sum_u coeff1_u S_u(x, var2)) * (sum_v coeff2_v S_v(x, var3))``
by summing ``schubmult_double(perm_dict1, v, var2, var3)`` scaled by
``coeff2_v`` over ``v`` in ``perm_dict2``.

<a id="schubmult.mult.double.schubmult_double"></a>

#### schubmult\_double

```python
def schubmult_double(perm_dict, v, var2=None, var3=None)
```

Multiply ``sum_u coeff_u S_u(x, var2)`` by the double Schubert polynomial ``S_v(x, var3)``.

Dispatches to the compiled ``schubmult_cpp`` kernel when available (and both
secondary alphabets are given), falling back to the pure-Python
implementation ``_schubmult_double_python`` otherwise.

**Arguments**:

- `perm_dict` - Mapping ``{Permutation: coeff}``.
- `v` - Permutation (or array-form list) indexing the Schubert polynomial to
  multiply by.
- `var2` - Secondary alphabet attached to ``perm_dict``'s permutations.
- `var3` - Secondary alphabet attached to ``v``.
  

**Returns**:

- `dict` - Coefficient dict ``{Permutation: coeff}`` (polynomials in ``var2``/``var3``).

<a id="schubmult.mult.double.schubmult_double_alt"></a>

#### schubmult\_double\_alt

```python
def schubmult_double_alt(perm_dict, v, var2=None, var3=None, index=1)
```

Alternate double Schubert product, built by peeling one variable of ``~v`` at a
time via ``pull_out_var`` instead of the ``theta``/v-path recursion.

Multiplies ``sum_u coeff_u S_u(x, var2)`` by ``S_v(x, var3)``, recursing on
``~new_v`` with the elementary symmetric factor coming from
``elem_sym_positional_perms`` at each step.

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_forwards"></a>

#### schubmult\_double\_alt\_from\_elems\_forwards

```python
def schubmult_double_alt_from_elems_forwards(perm_dict,
                                             v,
                                             var2=None,
                                             var3=None,
                                             index=1,
                                             elem_func=None)
```

``schubmult_double_alt`` generalized to an arbitrary elementary-symmetric-like
``elem_func(p, k, x_vars, y_vars)``, processing variables of ``~v`` from the first
pulled-out index forward.

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_backwards"></a>

#### schubmult\_double\_alt\_from\_elems\_backwards

```python
def schubmult_double_alt_from_elems_backwards(perm_dict,
                                              v,
                                              var2=None,
                                              var3=None,
                                              elem_func=None)
```

Like ``schubmult_double_alt_from_elems_forwards`` but processing ``~v``'s pulled-out
variables from the last descent backward, multiplying the elementary-symmetric
factor in *before* recursing (dispatches to the compiled kernel when available).

<a id="schubmult.mult.double.schubmult_double_alt_from_elems_backwards_backwards"></a>

#### schubmult\_double\_alt\_from\_elems\_backwards\_backwards

```python
def schubmult_double_alt_from_elems_backwards_backwards(
        perm_dict, v, var2=None, var3=None, elem_func=None)
```

Variant of ``_schubmult_double_alt_from_elems_backwards_python`` without the
per-``new_v`` memoization cache, recursing on the interim dict instead of the
original ``perm_dict`` at each pulled-out variable.

<a id="schubmult.mult.double.schubmult_double_from_elems"></a>

#### schubmult\_double\_from\_elems

```python
def schubmult_double_from_elems(perm_dict,
                                v,
                                var2=None,
                                var3=None,
                                elem_func=None)
```

``schubmult_double`` generalized to an arbitrary elementary-symmetric-like
``elem_func``, via the ``theta``/v-path recursion (rather than ``pull_out_var``).

Dispatches to the compiled kernel when available, falling back to
``_schubmult_double_from_elems_python``.

<a id="schubmult.mult.double.schubmult_double_down"></a>

#### schubmult\_double\_down

```python
def schubmult_double_down(perm_dict, v, var2=None, var3=None)
```

Down (descent) variant of ``_schubmult_double_python``, using ``elem_sym_perms_op``.

<a id="schubmult.mult.double.schub_coprod_double"></a>

#### schub\_coprod\_double

```python
def schub_coprod_double(mperm, indices, var2=None, var3=None)
```

Coproduct of the double Schubert polynomial ``S_mperm`` restricted to the
variable split named by ``indices``.

Analogue of ``schub_coprod_py``: multiplies the Grassmannian permutation for
``indices`` against ``mperm`` (via ``schubmult_double`` with a merged ``2N``
variable alphabet), splits each resulting permutation's window, and
substitutes the merged alphabet back to ``var2``/``var3``.

**Arguments**:

- `mperm` - Permutation (or array-form list) to take the coproduct of.
- `indices` - Iterable of 1-indexed positions selecting the variable split.
- `var2` - Secondary alphabet for the first factor's variables.
- `var3` - Secondary alphabet for the second factor's variables.
  

**Returns**:

- `dict` - Mapping ``{(firstperm, secondperm): coeff}``.

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

<a id="schubmult.utils.schub_lib"></a>

# schubmult.utils.schub\_lib

<a id="schubmult.utils.schub_lib.elem_sym_chains_groth"></a>

#### elem\_sym\_chains\_groth

```python
def elem_sym_chains_groth(orig_perm, p, k)
```

1 force mark, -1 force unmark, 0 otherwise

<a id="schubmult.utils.schub_lib.all_grassmannian_rc_graphs"></a>

#### all\_grassmannian\_rc\_graphs

```python
def all_grassmannian_rc_graphs(n: int, max_inv: int)
```

All RC graphs for Grassmannian permutations generated from partitions.

<a id="schubmult.utils.schub_lib.grassmannian_perms_from_partitions"></a>

#### grassmannian\_perms\_from\_partitions

```python
def grassmannian_perms_from_partitions(n: int,
                                       max_inv: int) -> list[Permutation]
```

Build Grassmannian permutations from partitions: pad to length n, reverse, then uncode.

<a id="schubmult.utils.schub_lib.partitions_with_sum_at_most"></a>

#### partitions\_with\_sum\_at\_most

```python
def partitions_with_sum_at_most(max_sum: int,
                                max_parts: int) -> list[tuple[int, ...]]
```

All partitions (weakly decreasing tuples) with at most max_parts parts and total <= max_sum.

<a id="schubmult.utils.schub_lib.hw_elementary_tensors"></a>

#### hw\_elementary\_tensors

```python
def hw_elementary_tensors(c, bounds=None)
```

Enumerate all sequences (A_1, ..., A_n) with A_i subset of [a_i] and |A_i| = c[i-1],
such that for every j the column sum d_j = #{i : j in A_i} is weakly decreasing in j.

Parameters
----------
c : sequence of nonnegative ints of length n
    A weak composition. Each c[i-1] must satisfy c[i-1] <= a_i.
bounds : sequence of positive ints of length n, optional
    A weakly increasing sequence a_1 <= a_2 <= ... <= a_n of positive integers.
    Defaults to (1, 2, ..., n).

Yields
------
tuple of tuples
    One tuple (A_1, ..., A_n) per valid family; each A_i is a sorted tuple of
    ints in [1, a_i] of length c[i-1].

<a id="schubmult.utils.tuple_utils"></a>

# schubmult.utils.tuple\_utils

<a id="schubmult.utils.tuple_utils.pad_tuple"></a>

#### pad\_tuple

```python
def pad_tuple(tup, length)
```

Pad a tuple-like with trailing zeros up to ``length``.

<a id="schubmult._scripts.count_sanity"></a>

# schubmult.\_scripts.count\_sanity

Count subwords of long_word(n) under two filters and compare.

<a id="schubmult._scripts.count_sanity.count_sylvester_subwords"></a>

#### count\_sylvester\_subwords

```python
def count_sylvester_subwords(code, n=None)
```

Count `1`: subwords of full (barred+unbarred) long word in Syl(F).

<a id="schubmult._scripts.count_sanity.count_reduced_with_forest"></a>

#### count\_reduced\_with\_forest

```python
def count_reduced_with_forest(code, n=None)
```

Count `2`: subwords of the FULL long word (barred+unbarred) whose value
sequence is a reduced word for perm = uncode(code) and whose
omega_invariant[0] (of reversed values) matches the principal-RC target.

<a id="schubmult._scripts.forest_graph_lr_rule"></a>

# schubmult.\_scripts.forest\_graph\_lr\_rule

<a id="schubmult._scripts.forest_graph_lr_rule.clear_rcgraph_caches"></a>

#### clear\_rcgraph\_caches

```python
def clear_rcgraph_caches()
```

Clear all known RCGraph and related caches to reduce memory overhead.

<a id="schubmult._scripts.quasi_dd_test"></a>

# schubmult.\_scripts.quasi\_dd\_test

<a id="schubmult._scripts.quasi_dd_test.t_index_iA"></a>

#### t\_index\_iA

```python
def t_index_iA(i, A)
```

Return j such that t_{i, A} = t_j, i.e. j = (\bar A)_i (1-indexed i-th elt of N \ A).

<a id="schubmult._scripts.quasi_dd_test.ev_A"></a>

#### ev\_A

```python
def ev_A(pol, x_arr, A, t_gen, n_xvars)
```

ev_A f = f(t_{(bar A)_1}, ..., t_{(bar A)_{n_xvars}}; t) — substitute x_j -> t_{(bar A)_j}.

<a id="schubmult._scripts.quasi_dd_test.factorization_from_code"></a>

#### factorization\_from\_code

```python
def factorization_from_code(code)
```

Columns factorization of forest F=1^{c_1}·2^{c_2}·... in Thompson monoid.

<a id="schubmult._scripts.quasi_dd_test.a_F_polynomial"></a>

#### a\_F\_polynomial

```python
def a_F_polynomial(f_poly, x_arr, t_gen, code, n_xvars)
```

Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
polynomial divided differences.

For factorization (i_1, ..., i_k) of F:
    a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
where A_j = i_{j+1} ⋆ ... ⋆ i_k.

<a id="schubmult._scripts.quasi_dd_test.schub_a_F_polynomial"></a>

#### schub\_a\_F\_polynomial

```python
def schub_a_F_polynomial(schub, x_arr, t_gen, code, n_xvars)
```

Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
polynomial divided differences.

For factorization (i_1, ..., i_k) of F:
    a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
where A_j = i_{j+1} ⋆ ... ⋆ i_k.

<a id="schubmult._scripts.quasi_dd_test.lrcoeff_a_F_polynomial"></a>

#### lrcoeff\_a\_F\_polynomial

```python
def lrcoeff_a_F_polynomial(schub, x_arr, t_gen, code, n_xvars)
```

Compute a_F(t) = [ev ⋆ E_F] f_poly  for forest with given code via direct
polynomial divided differences.

For factorization (i_1, ..., i_k) of F:
    a_F = ev_{A_0} ∘ E_{i_1, A_1} ∘ ... ∘ E_{i_k, ∅}  f
where A_j = i_{j+1} ⋆ ... ⋆ i_k.

<a id="schubmult._scripts.quasi_dd_test.enum_forest_codes"></a>

#### enum\_forest\_codes

```python
def enum_forest_codes(length, max_sum)
```

All forest codes (c_1, ..., c_length) with sum <= max_sum, c_i >= 0.

<a id="schubmult._scripts.graph_lr_rule"></a>

# schubmult.\_scripts.graph\_lr\_rule

<a id="schubmult._scripts.graph_lr_rule.clear_rcgraph_caches"></a>

#### clear\_rcgraph\_caches

```python
def clear_rcgraph_caches()
```

Clear all known RCGraph and related caches to reduce memory overhead.

<a id="schubmult._scripts.glide_lr_rule_nope"></a>

# schubmult.\_scripts.glide\_lr\_rule\_nope

<a id="schubmult._scripts.glide_lr_rule_nope.glide_rc_try"></a>

#### glide\_rc\_try

```python
def glide_rc_try(comp, beta, length)
```

GlidePoly polynomial of ``comp`` built by inverting the omega insertion.

We enumerate the compatible set-valued labelings ``kappa`` of the indexed
forest ``F = weak_composition_to_indfor(comp)`` (the glide definition), and
realize each labeling as a WCGraph by writing down its (word, compatible
sequence) pair *explicitly* -- the published inverse of the set-valued omega
insertion -- and feeding it to ``WCGraph.from_word_compatible``.

Concretely:

* The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
  labeling ``Q`` of the principal reduced RC graph form the omega pair for
  ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
  reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
  so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
  word position ``len(W) - Q(v)`` and carries the reduced letter
  ``ell(v) = W[len(W) - Q(v)]``.
* A labeling ``kappa`` places the letter ``ell(v)`` into every row
  ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
  ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
  a word whose rows are strictly decreasing, i.e. exactly the data
  ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
  values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
  equals ``x^kappa``.

Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
the degree ``-1`` homogenizer recording extra labels beyond one per node.

<a id="schubmult._scripts.grove_lr_rule_nope"></a>

# schubmult.\_scripts.grove\_lr\_rule\_nope

<a id="schubmult._scripts.grove_lr_rule_nope.grove_rc_try"></a>

#### grove\_rc\_try

```python
def grove_rc_try(comp, beta, length)
```

GrovePoly polynomial of ``comp`` built by inverting the omega insertion.

We enumerate the compatible set-valued labelings ``kappa`` of the indexed
forest ``F = weak_composition_to_indfor(comp)`` (the grove definition), and
realize each labeling as a WCGraph by writing down its (word, compatible
sequence) pair *explicitly* -- the published inverse of the set-valued omega
insertion -- and feeding it to ``WCGraph.from_word_compatible``.

Concretely:

* The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
  labeling ``Q`` of the principal reduced RC graph form the omega pair for
  ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
  reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
  so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
  word position ``len(W) - Q(v)`` and carries the reduced letter
  ``ell(v) = W[len(W) - Q(v)]``.
* A labeling ``kappa`` places the letter ``ell(v)`` into every row
  ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
  ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
  a word whose rows are strictly decreasing, i.e. exactly the data
  ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
  values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
  equals ``x^kappa``.

Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
the degree ``-1`` homogenizer recording extra labels beyond one per node.

<a id="schubmult._scripts.cauchy_dual_poly"></a>

# schubmult.\_scripts.cauchy\_dual\_poly

<a id="schubmult._scripts.cauchy_dual_poly.cauchy_dual_key"></a>

#### cauchy\_dual\_key

```python
def cauchy_dual_key(n)
```

Demazure atoms?

<a id="schubmult._scripts.cauchy_dual_poly.cauchy_dual_forest"></a>

#### cauchy\_dual\_forest

```python
def cauchy_dual_forest(n)
```

???

<a id="schubmult._scripts.check_omega_last_row"></a>

# schubmult.\_scripts.check\_omega\_last\_row

<a id="schubmult._scripts.check_omega_last_row.check_disjoint_vs_squash_omega"></a>

#### check\_disjoint\_vs\_squash\_omega

```python
def check_disjoint_vs_squash_omega(n: int) -> tuple[bool, int, int]
```

Compare omega_invariant[1] for disjoint_union vs squash_product.

Loop over minimal-length RC graphs rc. For each rc, loop over elementary
symmetric RC graphs of matching length and last descent, and verify:

    rc.disjoint_union(elem_sym_rc).omega_invariant[1]
    ==
    rc.squash_product(elem_sym_rc).omega_invariant[1]

<a id="schubmult._scripts.compu_double_forest"></a>

# schubmult.\_scripts.compu\_double\_forest

<a id="schubmult._scripts.compu_double_forest.class_key"></a>

#### class\_key

```python
def class_key(values)
```

Omega-class invariant: P-symbol of omega_insertion.

<a id="schubmult._scripts.compu_double_forest.build_match"></a>

#### build\_match

```python
def build_match(comp, n)
```

Match B-subwords to A-subwords by omega_insertion Q-symbol.

A := Sylvester subwords of long_word (canonical_forest == sylv_canon).
B := reduced-word subwords for uncode(comp) in the omega-class of the
     canonical reduced word of perm.

For each subword W (A or B), compute (P,Q) = omega_insertion(reverse(W)).
Both A and B subwords have shape F. Two subwords with identical Q place
their i-th reversed letter at the same node, so a_rev[i] <-> b_rev[i]
pairs them position-by-position; equivalently a[pos] <-> b[pos] in
original (long_word) order. The Q-bijection is unique and total.

<a id="schubmult._scripts.lr_rule_verify"></a>

# schubmult.\_scripts.lr\_rule\_verify

<a id="schubmult._scripts.lr_rule_verify.json_key_rep"></a>

#### json\_key\_rep

```python
def json_key_rep(k)
```

Serialize an arbitrary Python object `k` to a JSON-key-safe string.
Uses pickle + base64 with a 'PICKLE:' prefix; falls back to repr.

<a id="schubmult._scripts.lr_rule_verify.json_key_load"></a>

#### json\_key\_load

```python
def json_key_load(s)
```

Inverse of json_key_rep.

<a id="schubmult._scripts.lr_rule_verify.safe_save_recording"></a>

#### safe\_save\_recording

```python
def safe_save_recording(obj, filename, meta=None)
```

Save {meta: {...}, records: {json_key: value}} atomically.

<a id="schubmult._scripts.lr_rule_verify.safe_load_recording"></a>

#### safe\_load\_recording

```python
def safe_load_recording(filename)
```

Load the verification JSON.
Returns (records_dict, meta_dict).
Backwards compatible with old flat mapping format (treats file as records and meta={}).

<a id="schubmult._scripts.lr_rule_verify.recording_saver"></a>

#### recording\_saver

```python
def recording_saver(shared_recording_dict,
                    lock,
                    verification_filename,
                    stop_event,
                    meta=None,
                    sleep_time=10)
```

Background process that periodically snapshots `shared_recording_dict` and
writes it to `verification_filename`. The lock is held only briefly to
snapshot keys; the potentially slow file I/O and proxy lookups are done
outside the lock so workers are not blocked.

<a id="schubmult._scripts.lr_rule_verify.queue_producer"></a>

#### queue\_producer

```python
def queue_producer(task_queue, perms, n, num_processors, skip_id, irreducible,
                   base_level, sep_descs, dominant_only, w0_only,
                   shared_recording_dict, molevsagan)
```

Producer with periodic garbage collection.

<a id="schubmult._scripts.forest_lr_rule"></a>

# schubmult.\_scripts.forest\_lr\_rule

<a id="schubmult._scripts.forest_lr_rule.clear_rcgraph_caches"></a>

#### clear\_rcgraph\_caches

```python
def clear_rcgraph_caches()
```

Clear all known RCGraph and related caches to reduce memory overhead.

<a id="schubmult._scripts.elem_monom_formula"></a>

# schubmult.\_scripts.elem\_monom\_formula

<a id="schubmult._scripts.elem_monom_formula.elem_monom_formula"></a>

#### elem\_monom\_formula

```python
def elem_monom_formula(n)
```

???

<a id="schubmult._scripts.slide_lr_rule"></a>

# schubmult.\_scripts.slide\_lr\_rule

<a id="schubmult._scripts.slide_lr_rule.clear_rcgraph_caches"></a>

#### clear\_rcgraph\_caches

```python
def clear_rcgraph_caches()
```

Clear all known RCGraph and related caches to reduce memory overhead.

<a id="schubmult._scripts.pieri_formula_forest"></a>

# schubmult.\_scripts.pieri\_formula\_forest

<a id="schubmult._scripts.pieri_formula_forest.pieri_forest"></a>

#### pieri\_forest

```python
def pieri_forest(n)
```

???

<a id="schubmult._scripts.double_forest_polynomial"></a>

# schubmult.\_scripts.double\_forest\_polynomial

Compute double forest polynomials P_F(x; t) via the vine subword model.

Reference
---------
N. Bergeron, L. Gagnon, P. Nadeau, H. Spink, V. Tewari,
"Equivariant quasisymmetry and noncrossing partitions" (arXiv:2504.15234),
Theorem 5.1 and Section 5.1.

Given a code c = (c_1, c_2, ...) for an indexed forest F (so that c(F) = c),
and two sequences of symbols x = (x_1, x_2, ...) and t = (t_1, t_2, ...),
this returns the double forest polynomial

    P_F(x; t) = sum_{pi in R(omega_tilde_[n]; F)} wt(pi),

where omega_tilde_[n] is the long vine word and the sum is over its
subwords whose value-sequence is a Sylvester word for F.

<a id="schubmult._scripts.double_forest_polynomial.forest_from_code"></a>

#### forest\_from\_code

```python
def forest_from_code(code)
```

Build the indexed forest F with c(F) = ``code`` via the Thompson
monoid factorization F = 1^{c_1} . 2^{c_2} . 3^{c_3} . ...

<a id="schubmult._scripts.double_forest_polynomial.sylvester_word"></a>

#### sylvester\_word

```python
def sylvester_word(forest)
```

Return one Sylvester word of ``forest`` via pre-order traversal
(root, left subtree, right subtree) over each non-trivial tree.

<a id="schubmult._scripts.double_forest_polynomial.long_word"></a>

#### long\_word

```python
def long_word(n)
```

Build omega_tilde_[n] as a list of letter records.

omega^(k)_[n] = (n, n-1, ..., k+1, k, k+1_bar, ..., n-1_bar, n_bar)
Concatenate for k = 1, ..., n.

<a id="schubmult._scripts.double_forest_polynomial.letter_weight"></a>

#### letter\_weight

```python
def letter_weight(letter, x_gen, t_gen)
```

wt(j^(i))      = x_i - t_j   (unbarred)
wt(j_bar^(i)) = t_j - t_i   (barred)

<a id="schubmult._scripts.double_forest_polynomial.double_forest_polynomial"></a>

#### double\_forest\_polynomial

```python
def double_forest_polynomial(code, x_gen, t_gen, n=None)
```

Compute P_F(x; t) for the indexed forest F with c(F) = ``code``.

Parameters
----------
code : sequence[int]
    The code (c_1, c_2, ...) of F. Trailing zeros are ignored.
x_gen : callable[int -> Expr]
    Maps i -> x_i  (the non-equivariant variables).
t_gen : callable[int -> Expr]
    Maps j -> t_j  (the equivariant variables).
n : int, optional
    Use the long word omega_tilde_[n]. If ``None``, the smallest n with
    F supported in {1, ..., n+1} is used.

Returns
-------
sympy.Expr
    The expanded polynomial P_F(x; t).

<a id="schubmult._scripts.double_forest_polynomial.decompose_double_forest_tensor"></a>

#### decompose\_double\_forest\_tensor

```python
def decompose_double_forest_tensor(poly, x_genset, t_genset, length)
```

Decompose a double forest polynomial into t-forest ⊗ x-forest terms.

Returns a dict mapping ``(t_code, x_code) -> integer/symbolic coefficient``.

<a id="schubmult._scripts.double_forest_polynomial.DoubleForestPolynomialBasis"></a>

## DoubleForestPolynomialBasis Objects

```python
class DoubleForestPolynomialBasis()
```

Concrete basis for double forest polynomials.

A basis element is indexed by a pair ``(t_code, x_code)`` and represents
``Forest_t(t_code) ⊗ Forest_x(x_code)``.

<a id="schubmult._scripts.double_forest_polynomial.ForestDoubleElement"></a>

## ForestDoubleElement Objects

```python
class ForestDoubleElement(dict)
```

Formal element in the abstract double-forest basis.

Stored as {forest_code: coeff}, representing
    sum_F coeff[F] * DF(F),
where DF(F) is an abstract basis symbol (not expanded into x/t variables).

<a id="schubmult._scripts.double_forest_polynomial.ForestDoubleElement.expand_polynomial"></a>

#### expand\_polynomial

```python
def expand_polynomial()
```

Expand abstract basis expression to a concrete x/t polynomial.

<a id="schubmult._scripts.double_forest_polynomial.ForestDouble"></a>

## ForestDouble Objects

```python
class ForestDouble()
```

Abstract double-forest basis algebra.

Basis symbols are indexed by one forest code F and represent P_F(x; t)
abstractly. Coefficients live in the polynomial ring of equivariant vars.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_coefficients"></a>

#### extract\_double\_forest\_coefficients

```python
def extract_double_forest_coefficients(poly, x_genset, length)
```

Return all coefficients a_F(t) in f = sum_F a_F(t) P_F(x; t).

Implemented by iterative top-degree peeling:
1) convert the current residual into ForestPolyBasis in x,
2) read top-degree coefficients,
3) subtract those coefficients times the corresponding double forest basis
   polynomials P_F,
4) repeat until the residual vanishes.

This avoids treating the double-forest basis as a plain one-shot x-basis
conversion and follows the triangular subtraction strategy.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_coefficient"></a>

#### extract\_double\_forest\_coefficient

```python
def extract_double_forest_coefficient(poly, forest_code, x_genset, length)
```

Return the single coefficient a_F(t) of P_F in f.

Equivalent to the operator formula a_F(t) = [ev star e_F] f from
Theorem 10.9 (Coefficient extraction and star-composition).

Note: this returns a t-polynomial coefficient in the x-forest basis.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_tensor_coefficients"></a>

#### extract\_double\_forest\_tensor\_coefficients

```python
def extract_double_forest_tensor_coefficients(poly, x_genset, t_genset,
                                              length)
```

Return true double-forest tensor-basis coefficients.

Expands poly as
    poly = sum_{A,B} c_{A,B} Forest_t(A) \otimes Forest_x(B),
returning a dict mapping (A, B) -> c_{A,B}.

<a id="schubmult._scripts.double_forest_polynomial.extract_double_forest_tensor_coefficient"></a>

#### extract\_double\_forest\_tensor\_coefficient

```python
def extract_double_forest_tensor_coefficient(poly, t_code, x_code, x_genset,
                                             t_genset, length)
```

Return c_{A,B} for one tensor basis pair A=t_code, B=x_code.

<a id="schubmult._scripts.double_forest_polynomial.monk_style_degree_one_product"></a>

#### monk\_style\_degree\_one\_product

```python
def monk_style_degree_one_product(code,
                                  i,
                                  x_gen,
                                  t_gen,
                                  x_genset,
                                  t_genset,
                                  n=None,
                                  length=None)
```

Compute a Monk-style decomposition for degree-1 double forest multiplication.

The product
    P_{e_i}(x; t) * P_code(x; t)
is decomposed as
    diagonal_coeff(t) * P_code(x; t) + sum_{beta != code} c_beta(t) * P_beta(x; t),
and also expanded in the tensor basis t-forest ⊗ x-forest.

<a id="schubmult._scripts.grothmult_double"></a>

# schubmult.\_scripts.grothmult\_double

<a id="schubmult._scripts.grothmult_double.groth_posify"></a>

#### groth\_posify

```python
def groth_posify(val, var2, var3, msg)
```

Positive FGL form of a coefficient, at ``beta = 1``.

Basis: all products of ``z_i - y_j``, ``(1 + y_j)^{-1}``, ``(1 + z_i)^{-1}``.
Substituting ``u_j = 1 + y_j``, ``v_i = 1 + z_i`` (so ``y = u - 1``,
``z = v - 1``) turns the value into a Laurent polynomial in ``u, v``; clearing
the monomial denominator makes everything polynomial.  Candidates are
difference products times leftover clearing monomials, expanded once, with
``as_coefficients_dict`` terms used as opaque basis vectors for an integer
LP.  ``beta`` is restored as ``beta**(`diffs` - d)`` with the Laurent atoms
``(1 + beta*y)`` degree 0 by construction.

<a id="schubmult.symbolic.poly.schub_poly"></a>

# schubmult.symbolic.poly.schub\_poly

<a id="schubmult.symbolic.poly.schub_poly.isobaric_strip_on_dschub_dict"></a>

#### isobaric\_strip\_on\_dschub\_dict

```python
def isobaric_strip_on_dschub_dict(start, length, perm_dict, coeff_genset,
                                  beta)
```

Apply one isobaric strip to a whole ``{perm: coeff}`` dict, folded.

Coefficients landing on the same permutation merge at every stage instead of
being carried per input basis element, mirroring ``compute_vpathdicts``.

<a id="schubmult.symbolic.poly.schub_poly.apply_isobaric_to_schub_dict"></a>

#### apply\_isobaric\_to\_schub\_dict

```python
def apply_isobaric_to_schub_dict(diff_perm, perm_dict, coeff_genset, beta)
```

Fold every strip of ``diff_perm`` over the whole dict, merging between strips.

<a id="schubmult.symbolic.poly.schub_poly.groth_elem_as_schub_dict"></a>

#### groth\_elem\_as\_schub\_dict

```python
@cache
def groth_elem_as_schub_dict(perm, beta)
```

Expand a Grothendieck basis element G_perm into Schubert basis terms.

Uses the strip-isobaric construction from the groth_elem_as_schub method.
Returns ``{Permutation: coeff}``.

<a id="schubmult.symbolic.poly.variables"></a>

# schubmult.symbolic.poly.variables

<a id="schubmult.symbolic.poly.variables.GeneratingSet"></a>

## GeneratingSet Objects

```python
class GeneratingSet(GeneratingSet_base)
```

<a id="schubmult.symbolic.poly.variables.GeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet"></a>

## MaskedGeneratingSet Objects

```python
class MaskedGeneratingSet(GeneratingSet_base)
```

<a id="schubmult.symbolic.poly.variables.MaskedGeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.CustomGeneratingSet"></a>

## CustomGeneratingSet Objects

```python
class CustomGeneratingSet(GeneratingSet_base)
```

<a id="schubmult.symbolic.poly.variables.CustomGeneratingSet.__call__"></a>

#### \_\_call\_\_

```python
def __call__(index)
```

1-indexed

<a id="schubmult.symbolic.poly.variables.genset_dict_from_expr"></a>

#### genset\_dict\_from\_expr

```python
def genset_dict_from_expr(expr, genset, length=None)
```

Transform expressions into a multinomial form given generators.

