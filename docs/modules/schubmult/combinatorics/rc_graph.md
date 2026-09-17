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

