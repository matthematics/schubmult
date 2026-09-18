<a id="schubmult.rings.quasisymmetric_functions"></a>

# schubmult.rings.quasisymmetric\_functions

`QSym`: quasisymmetric functions in the monomial basis ``M_alpha``.

Keys are compositions; the product is the quasi-shuffle (stuffle) of compositions, and
``expand(n)`` gives the monomial quasisymmetric polynomial in ``n`` variables. `QSym.quasi_schur`
builds quasi-Schur functions by enumerating standard composition tableaux.

<a id="schubmult.rings.quasisymmetric_functions.monomial_quasisym"></a>

#### monomial\_quasisym

```python
def monomial_quasisym(comp, length, genset)
```

The monomial quasisymmetric polynomial ``M_comp(x_1, ..., x_length)``: the sum of
``x_{i_1}^{c_1} ... x_{i_k}^{c_k}`` over ``i_1 < ... < i_k <= length``, built by recursion on
whether ``x_length`` is used. Zero if ``comp`` contains a zero part.

<a id="schubmult.rings.quasisymmetric_functions.stuffle"></a>

#### stuffle

```python
def stuffle(alpha, beta)
```

The quasi-shuffle (stuffle) product of two compositions: at each step take the first part
of ``alpha``, the first part of ``beta``, or their sum. Returns ``{composition: coeff}``;
this is the product rule of the monomial basis ``M_alpha M_beta``.

<a id="schubmult.rings.quasisymmetric_functions.quasi_schur_to_monomial"></a>

#### quasi\_schur\_to\_monomial

```python
def quasi_schur_to_monomial(comp)
```

Monomial-basis expansion of the quasi-Schur function of shape ``comp``: counts standard
composition tableaux (rows strictly increasing, columns weakly increasing) of that shape by
the descent composition of their row reading word. Enumerates all ``n!`` fillings, so only
small shapes are practical.

<a id="schubmult.rings.quasisymmetric_functions.QSymElement"></a>

## QSymElement Objects

```python
class QSymElement(BaseSchubertElement)
```

Element of `QSym`: a dict from compositions to coefficients in the monomial basis.

<a id="schubmult.rings.quasisymmetric_functions.QSymElement.expand"></a>

#### expand

```python
def expand(num_vars)
```

The quasisymmetric polynomial in ``num_vars`` variables of the ring's generating set.

<a id="schubmult.rings.quasisymmetric_functions.QSym"></a>

## QSym Objects

```python
class QSym(BaseSchubertRing)
```

Quasisymmetric functions in the monomial basis; ``QSym()(2, 1)`` is ``M_(2,1)``. See the module docstring.

<a id="schubmult.rings.quasisymmetric_functions.QSym.mul_pair"></a>

#### mul\_pair

```python
def mul_pair(a, b)
```

Product of two basis compositions: the `stuffle`.

<a id="schubmult.rings.quasisymmetric_functions.QSym.mul"></a>

#### mul

```python
def mul(a, b)
```

Bilinear extension of `mul_pair`.

<a id="schubmult.rings.quasisymmetric_functions.QSym.printing_term"></a>

#### printing\_term

```python
def printing_term(comp)
```

Display as ``Mx(alpha)`` (label from the generating set).

<a id="schubmult.rings.quasisymmetric_functions.QSym.new"></a>

#### new

```python
def new(*x)
```

The basis element ``M_x`` for the composition given as positional parts.

<a id="schubmult.rings.quasisymmetric_functions.QSym.quasi_schur"></a>

#### quasi\_schur

```python
def quasi_schur(*comp)
```

The quasi-Schur function of shape ``comp`` in the monomial basis (see `quasi_schur_to_monomial`).

>>> QS = QSym()
>>> QS.quasi_schur(2, 1)

