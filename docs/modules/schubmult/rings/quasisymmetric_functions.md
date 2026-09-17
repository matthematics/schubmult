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

