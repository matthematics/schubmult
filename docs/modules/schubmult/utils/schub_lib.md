<a id="schubmult.utils.schub_lib"></a>

# schubmult.utils.schub\_lib

Combinatorial kernels behind the Schubert multiplication algorithms.

The central objects are the *v-path dictionaries* (`compute_vpathdicts`): for a theta code
``th`` and target permutation ``vmu`` they record, level by level, the Bruhat-descent moves
(`kdown_perms`) along which the Schubert product is accumulated in `schubmult.mult.single` and
its double/quantum relatives. Also here: the ``elem_sym_perms*``/``complete_sym_perms*`` families
(Pieri-type moves for multiplying by elementary/complete symmetric polynomials, including quantum
and Grothendieck variants), `pull_out_var` (the recursion behind `schubpoly`), coefficient
reductions (`reduce_coeff`, `reduce_descents`, ``try_reduce_*``) that shrink a triple ``(u, v, w)``
before computing an LR coefficient, and the parabolic ``q``-vector checks.

<a id="schubmult.utils.schub_lib.double_elem_sym_q"></a>

#### double\_elem\_sym\_q

```python
def double_elem_sym_q(u, p1, p2, k, q_var=q_var)
```

Pairs of consecutive quantum Pieri moves ``u -> perm1 -> perm2`` (degrees ``p1`` then ``p2`` in ``k``
variables) whose cycles are compatible; returned as ``{(perm1, udiff1, q1): [(perm2, udiff2, q2), ...]}``.

<a id="schubmult.utils.schub_lib.will_formula_work"></a>

#### will\_formula\_work

```python
def will_formula_work(u, v)
```

Whether the fast dominant-descent formula applies to the pair ``(u, v)``: ``v^{-1} * mu_v`` must
already be the identity (no descents to reduce).

<a id="schubmult.utils.schub_lib.try_reduce_u"></a>

#### try\_reduce\_u

```python
def try_reduce_u(u, v, w)
```

Move a zero in the code of ``u`` past a nonzero entry by swapping adjacent positions in ``(u, w)``
or ``(u, v)`` where the LR coefficient is unchanged, aiming for ``u`` to one-dominate ``w``.

<a id="schubmult.utils.schub_lib.reduce_descents"></a>

#### reduce\_descents

```python
def reduce_descents(u, v, w)
```

Cancel common descents of ``w`` with ``v`` (or ``u``) by simultaneous adjacent swaps, stopping once
one of the fast-path conditions holds.

<a id="schubmult.utils.schub_lib.is_reducible"></a>

#### is\_reducible

```python
def is_reducible(v)
```

Whether the code of ``v`` has no nonzero entry after a zero (so ``v`` is a shifted dominant-like shape).

<a id="schubmult.utils.schub_lib.try_reduce_v"></a>

#### try\_reduce\_v

```python
def try_reduce_v(u, v, w)
```

`try_reduce_u` with the roles of ``u`` and ``v`` exchanged, aiming for `is_reducible`.

<a id="schubmult.utils.schub_lib.reduce_coeff"></a>

#### reduce\_coeff

```python
def reduce_coeff(u, v, w)
```

Reduce the LR triple ``(u, v, w)`` by dividing out the dominant permutations of the theta codes of
``u^{-1}`` and ``v^{-1}``: if ``w`` decomposes compatibly, return the smaller triple
``(u mu_A, v mu_B, w')`` with the same coefficient; otherwise return the input.

<a id="schubmult.utils.schub_lib.pull_out_var"></a>

#### pull\_out\_var

```python
def pull_out_var(vnum, v)
```

All ways to factor the variable ``x_vnum`` out of ``S_v``: returns pairs ``(indices, v')`` such that
``S_v = sum prod_{p in indices} (x_vnum - y_p) * S_{v'}`` with ``x_vnum`` absent from ``S_{v'}``.

<a id="schubmult.utils.schub_lib.divdiffable"></a>

#### divdiffable

```python
def divdiffable(v, u)
```

``v u^{-1}`` if ``u <= v`` with lengths adding (so ``partial_{u}`` can be applied), else ``[]``.

<a id="schubmult.utils.schub_lib.kdown_perms"></a>

#### kdown\_perms

```python
def kdown_perms(perm: Permutation, monoperm: Permutation, p: int,
                k: int) -> list[tuple[Permutation, int, int]]
```

One level of the v-path recursion: all ``(new_perm, pp, sign)`` obtained from ``perm`` by up to ``p``
Bruhat descents through position ``k`` (swapping position ``k-1`` with an untouched position on
either side) such that ``new_perm * monoperm`` has the expected length.

<a id="schubmult.utils.schub_lib.rc_graph_set"></a>

#### rc\_graph\_set

```python
def rc_graph_set(perm)
```

All RC graphs of ``perm`` as ``(row_labels, reduced_word)`` pairs, built recursively with `pull_out_var`.

<a id="schubmult.utils.schub_lib.compute_vpathdicts_cached"></a>

#### compute\_vpathdicts\_cached

```python
@cache
def compute_vpathdicts_cached(th, vmu)
```

Build the v-path dictionaries for theta code ``th`` and target ``vmu``.

Working from the top level down, each level ``i`` maps a permutation to the `kdown_perms`
moves available at that level; the result is re-indexed as
``vpathdicts[i][source] = {(target, degree, sign), ...}`` for forward accumulation.

<a id="schubmult.utils.schub_lib.compute_vpathdicts"></a>

#### compute\_vpathdicts

```python
def compute_vpathdicts(th, vmu)
```

Cached `compute_vpathdicts_cached` accepting a list ``th``.

<a id="schubmult.utils.schub_lib.check_blocks"></a>

#### check\_blocks

```python
def check_blocks(qv, parabolic_index)
```

Parabolic admissibility of a ``q``-exponent vector: over every interval within each block of
consecutive parabolic indices, the sum of `omega` values must be 0 or -1.

<a id="schubmult.utils.schub_lib.reduce_q_coeff"></a>

#### reduce\_q\_coeff

```python
def reduce_q_coeff(u, v, w, qv)
```

One quantum reduction step: find a position where swapping adjacent entries of ``v`` (or ``u``)
together with ``w`` -- adjusting the ``q``-vector when ``w`` had no descent there -- leaves the
quantum LR coefficient unchanged. Returns ``(u, v, w, qv, changed)``.

<a id="schubmult.utils.schub_lib.reduce_q_coeff_u_only"></a>

#### reduce\_q\_coeff\_u\_only

```python
def reduce_q_coeff_u_only(u, v, w, qv)
```

`reduce_q_coeff` restricted to swaps on ``u``.

<a id="schubmult.utils.schub_lib.elem_sym_perms_q"></a>

#### elem\_sym\_perms\_q

```python
def elem_sym_perms_q(orig_perm, p, k, q_var=q_var)
```

Quantum Pieri moves for ``E_p^q(x_1..x_k)``: all ``(perm, degree, q_monomial)`` reachable from
``orig_perm`` by up to ``p`` transpositions ``(i, j)`` with ``i < k <= j`` on still-untouched positions
``i``; a Bruhat *descent* (cyclic quantum move) contributes ``q_{i+1} ... q_j``.

<a id="schubmult.utils.schub_lib.elem_sym_perms_q_op"></a>

#### elem\_sym\_perms\_q\_op

```python
def elem_sym_perms_q_op(orig_perm, p, k, n, q_var=q_var)
```

Adjoint (downward) version of `elem_sym_perms_q` on permutations padded to length ``n``.

<a id="schubmult.utils.schub_lib.elem_sym_perms"></a>

#### elem\_sym\_perms

```python
def elem_sym_perms(orig_perm, p, k)
```

Pieri moves for ``e_p(x_1..x_k)``: all ``(perm, degree)`` reachable from ``orig_perm`` by up to ``p``
Bruhat ascents ``(i, j)`` with ``i < k <= j``, each position ``i`` used at most once (values
swapped in strictly decreasing order to avoid double counting).

<a id="schubmult.utils.schub_lib.elem_sym_perms_groth"></a>

#### elem\_sym\_perms\_groth

```python
def elem_sym_perms_groth(orig_perm, p, k)
```

Grothendieck variant of `elem_sym_perms`: positions in ``i < k`` may be reused, with ``j``
weakly decreasing along the chain.

<a id="schubmult.utils.schub_lib.elem_sym_chains_groth"></a>

#### elem\_sym\_chains\_groth

```python
def elem_sym_chains_groth(orig_perm, p, k)
```

Marked Bruhat chains for the Grothendieck Pieri rule: chains of ascents through position ``k``
together with a marking vector (1 force mark, -1 force unmark, 0 free) enforcing the ordering
conditions on consecutive covers.

<a id="schubmult.utils.schub_lib.elem_sym_positional_perms"></a>

#### elem\_sym\_positional\_perms

```python
def elem_sym_positional_perms(orig_perm, p, *k)
```

Pieri moves for ``e_p`` in an arbitrary set of variable positions ``k`` (1-indexed): Bruhat ascents
swapping an untouched position in ``k`` with a position outside ``k``, returned as
``(perm, degree, sign)`` with sign ``-1`` when the outside position is to the left.

<a id="schubmult.utils.schub_lib.elem_sym_positional_perms_q"></a>

#### elem\_sym\_positional\_perms\_q

```python
def elem_sym_positional_perms_q(orig_perm, p, *k, q_var=q_var)
```

Quantum version of `elem_sym_positional_perms`: returns ``(perm, degree, sign, q_monomial)``.

<a id="schubmult.utils.schub_lib.complete_sym_perms_op"></a>

#### complete\_sym\_perms\_op

```python
def complete_sym_perms_op(orig_perm, p, k)
```

Downward Pieri moves for ``h_p(x_1..x_k)``: ``{perm: [(i, j), ...]}`` recording the chain of Bruhat
descents (position ``i < k`` with an untouched position ``j``).

<a id="schubmult.utils.schub_lib.complete_sym_perms"></a>

#### complete\_sym\_perms

```python
def complete_sym_perms(orig_perm, p, k)
```

Pieri moves for ``h_p(x_1..x_k)``: ``{perm: [(i, j), ...]}`` recording the chain of Bruhat ascents
(position ``i < k`` may repeat; ``j`` untouched).

<a id="schubmult.utils.schub_lib.complete_sym_positional_perms"></a>

#### complete\_sym\_positional\_perms

```python
def complete_sym_positional_perms(orig_perm, p, *k)
```

`complete_sym_perms` for an arbitrary set of positions ``k`` (1-indexed), with signs as in
`elem_sym_positional_perms`; returns ``(perm, degree, sign)`` triples.

<a id="schubmult.utils.schub_lib.elem_sym_perms_op"></a>

#### elem\_sym\_perms\_op

```python
def elem_sym_perms_op(orig_perm, p, k)
```

Downward (Bruhat descent) version of `elem_sym_perms`.

<a id="schubmult.utils.schub_lib.is_split_two"></a>

#### is\_split\_two

```python
def is_split_two(u, v, w)
```

For ``inv(w) - inv(u) == 2``: whether ``v^{-1} w`` is a product of two disjoint cycles; returns
``(True, cycles)`` or ``(False, [])``.

<a id="schubmult.utils.schub_lib.is_coeff_irreducible"></a>

#### is\_coeff\_irreducible

```python
def is_coeff_irreducible(u, v, w)
```

Whether none of the fast paths / reductions applies to the LR triple ``(u, v, w)``.

<a id="schubmult.utils.schub_lib.is_hook"></a>

#### is\_hook

```python
def is_hook(cd)
```

Whether the code ``cd`` is a hook shape: a run of 1's optionally followed by a single larger entry, then zeros.

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

<a id="schubmult.utils.schub_lib.fff"></a>

#### fff

```python
def fff(chain)
```

Number of forced marks (``1``) in a marked chain from `elem_sym_chains_groth`.

<a id="schubmult.utils.schub_lib.ppp"></a>

#### ppp

```python
def ppp(chain)
```

Number of forced unmarks (``-1``) in a marked chain from `elem_sym_chains_groth`.

<a id="schubmult.utils.schub_lib.groth_pieri_mul"></a>

#### groth\_pieri\_mul

```python
def groth_pieri_mul(perm_dict, p, kk, beta)
```

Grothendieck Pieri rule: multiply the Grothendieck expansion ``perm_dict`` by ``G`` of the
Grassmannian permutation of ``e_p(x_1..x_kk)``. Sums over marked chains from
`elem_sym_chains_groth` with weight ``beta^(length - p) * C(n, k)`` where the free marks are
chosen binomially.

