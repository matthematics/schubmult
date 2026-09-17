<a id="schubmult.combinatorics.chute_move_element"></a>

# schubmult.combinatorics.chute\_move\_element

Chute moves on RC graphs: a marked-row wrapper tracking a before/after pair of RC graphs
related by simultaneous chute moves on a chosen set of non-adjacent rows.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement"></a>

## ChuteMoveElement Objects

```python
class ChuteMoveElement(GridPrint)
```

An RC graph together with the result of applying chute moves at ``rows``.

``rows`` must be pairwise non-adjacent; for each row, moves the element hanging
off the row's end down into the first available gap in the row below,
raising ``ValueError`` if no valid chute move exists. Stores the pair
``(original, moved)`` RC graphs.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.product"></a>

#### product

```python
def product(other)
```

Stack ``self`` above ``other`` (via the underlying RC graph product) and combine their marked rows.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.chute_degree"></a>

#### chute\_degree

```python
@property
def chute_degree()
```

Number of marked rows (simultaneous chute moves applied).

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.chute_move_rows"></a>

#### chute\_move\_rows

```python
@property
def chute_move_rows()
```

The set of marked row indices.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.cols"></a>

#### cols

```python
@property
def cols()
```

Number of columns of the underlying RC graph.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.rows"></a>

#### rows

```python
@property
def rows()
```

Number of rows of the underlying RC graph.

<a id="schubmult.combinatorics.chute_move_element.ChuteMoveElement.print_element"></a>

#### print\_element

```python
@property
def print_element()
```

A `GridPrint`-compatible view highlighting the before/after cells at the marked rows.

