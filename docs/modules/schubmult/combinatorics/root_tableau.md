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

