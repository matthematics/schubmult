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

