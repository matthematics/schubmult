<a id="schubmult.combinatorics.pipe_dream"></a>

# schubmult.combinatorics.pipe\_dream

Classical pipe dreams (rectangular grid of crosses/bumps encoding a reduced word),
with conversion to/from `RCGraph`/`WCGraph` and grid symmetries (transpose-like duals,
inversion, vertical reflection).

<a id="schubmult.combinatorics.pipe_dream.PipeDream"></a>

## PipeDream Objects

```python
class PipeDream(PlanarHistory, GridPrint)
```

A pipe dream: pipes enter from the left and exit at the top of a triangular grid of
``CROSS``/``BUMP`` tiles. ``perm`` recovers the induced permutation; ``perm_word`` its
associated (not necessarily reduced) word.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.perm_word"></a>

#### perm\_word

```python
@property
def perm_word()
```

The word read off the crosses in native (left-to-right, top-to-bottom) orientation.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.is_reduced"></a>

#### is\_reduced

```python
@property
def is_reduced()
```

Whether ``perm_word`` is a reduced word for ``perm`` (length equals the number of crosses).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.perm"></a>

#### perm

```python
@property
def perm()
```

The permutation induced by this pipe dream (0-Hecke product of the cross word).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.to_rc_graph"></a>

#### to\_rc\_graph

```python
def to_rc_graph()
```

Convert to an `RCGraph` (crosses of row ``i`` become that row's column labels).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.to_wc_graph"></a>

#### to\_wc\_graph

```python
def to_wc_graph()
```

Convert to a `WCGraph`, analogous to ``to_rc_graph``.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.from_rc_graph"></a>

#### from\_rc\_graph

```python
@classmethod
def from_rc_graph(cls, rc_graph)
```

Build the pipe dream whose crosses are exactly the RC graph's marked positions.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.from_wc_graph"></a>

#### from\_wc\_graph

```python
@classmethod
def from_wc_graph(cls, wc_graph)
```

Build the pipe dream whose crosses are exactly the WC graph's marked positions.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.co_pipe_dream"></a>

#### co\_pipe\_dream

```python
def co_pipe_dream()
```

Swap crosses and bumps under the anti-diagonal reflection (the "co" dual pipe dream).

<a id="schubmult.combinatorics.pipe_dream.PipeDream.co_stinkbat_pipe_dream"></a>

#### co\_stinkbat\_pipe\_dream

```python
def co_stinkbat_pipe_dream()
```

Like ``co_pipe_dream`` but preserving (rather than swapping) cross/bump identity under the reflection.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.inverse_pipe_dream"></a>

#### inverse\_pipe\_dream

```python
def inverse_pipe_dream()
```

Pipe dream for ``~self.perm``, obtained by reflecting bumps/crosses through the anti-diagonal.

<a id="schubmult.combinatorics.pipe_dream.PipeDream.reflect_vertically"></a>

#### reflect\_vertically

```python
def reflect_vertically()
```

Flip the grid top-to-bottom.

