<a id="schubmult.combinatorics.set_word"></a>

# schubmult.combinatorics.set\_word

Set-valued crystal words: `SetLetter` (a subset-of-``{1..n}`` letter with a
GL_n-type crystal structure) and `SetWord` (a tensor of such letters), with
conversions to/from `WCGraph`.

<a id="schubmult.combinatorics.set_word.SetLetter"></a>

## SetLetter Objects

```python
class SetLetter(CrystalGraph, frozenset)
```

A set of ints with a sqrt(gl_n) crystal structure.

<a id="schubmult.combinatorics.set_word.SetLetter.crystal_length"></a>

#### crystal\_length

```python
def crystal_length()
```

The ambient rank ``n`` (number of crystal indices).

<a id="schubmult.combinatorics.set_word.SetLetter.crystal_weight"></a>

#### crystal\_weight

```python
@property
def crystal_weight()
```

Weight vector: multiplicity of each value ``1..n`` in the set.

<a id="schubmult.combinatorics.set_word.SetLetter.raising_operator"></a>

#### raising\_operator

```python
def raising_operator(i)
```

``e_i``: move an element from ``i+1`` to ``i`` if that increases the weight at ``i``, else ``None``.

<a id="schubmult.combinatorics.set_word.SetLetter.lowering_operator"></a>

#### lowering\_operator

```python
def lowering_operator(i)
```

``f_i``: the inverse move to ``raising_operator``, or ``None`` if undefined.

<a id="schubmult.combinatorics.set_word.SetWord"></a>

## SetWord Objects

```python
class SetWord(CrystalGraphTensor)
```

A tuple of SetLetters with a sqrt(gl_n) crystal structure.

<a id="schubmult.combinatorics.set_word.SetWord.to_wc_graph"></a>

#### to\_wc\_graph

```python
def to_wc_graph(rows)
```

Convert to a `WCGraph` with the given number of rows: column ``j`` gets a reflection
at each row in ``self.factors[j-1]``.

<a id="schubmult.combinatorics.set_word.SetWord.from_wc_graph"></a>

#### from\_wc\_graph

```python
@classmethod
def from_wc_graph(cls, wc)
```

Inverse of ``to_wc_graph``: build a `SetWord` from a `WCGraph`, one `SetLetter` per column.

