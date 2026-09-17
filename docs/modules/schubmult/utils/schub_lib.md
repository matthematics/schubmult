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

