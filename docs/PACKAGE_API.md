# schubmult — Python package API

Practical reference for using **schubmult** as a Python library. For the CLI
tools shipped with the package see the
[README](https://github.com/matthematics/schubmult).

This document targets schubmult **4.0.x**.

---

## Installation

```bash
pip install schubmult                # PyPI release (when available)
# or, for a development checkout:
git clone https://github.com/matthematics/schubmult
cd schubmult && pip install -e .
```


## Conceptual map

| Layer | What it is | Where it lives |
|---|---|---|
| **Combinatorics** | `Permutation`, RC-graphs, root tableaux | [`schubmult.combinatorics`](https://github.com/matthematics/schubmult/tree/main/src/schubmult/combinatorics) |
| **Multiplication kernels** | Pure functions `schubmult_*` on dicts | [`schubmult.mult`](https://github.com/matthematics/schubmult/tree/main/src/schubmult/mult) |
| **Rings** | Algebraic ring objects with `+`, `*`, `expand()` | [`schubmult.rings`](https://github.com/matthematics/schubmult/tree/main/src/schubmult/rings) |
| **Symbolic** | Hybrid SymPy/SymEngine facade and variable factories | [`schubmult.symbolic`](https://github.com/matthematics/schubmult/tree/main/src/schubmult/symbolic) |
| **CLI scripts** | `schubmult_py`, `schubmult_double`, `schubmult_q`, `schubmult_q_double` | [`schubmult._scripts`](https://github.com/matthematics/schubmult/tree/main/src/schubmult/_scripts) |

The two everyday entry points are the **rings** (`Sx`, `DSx`, …) for high-level
algebraic work and the **kernels** (`schubmult_py`, `schubmult_double`, …) for
fast dict-level computations.

---

## `Permutation` — the central object

```python
from schubmult import Permutation, uncode

w = Permutation([3, 1, 2])     # one-line notation: w(1)=3, w(2)=1, w(3)=2
w.code                         # Lehmer code (full)
w.trimcode                     # Lehmer code with trailing zeros stripped
w.descents()                   # descent set (0-indexed by default)
w.inv                          # number of inversions = length
w.inversion_set()              # set of inversion pairs (i,j)
w.is_dominant                  # bool
w.is_vexillary                 # bool
w.shape                        # partition shape (when dominant/vexillary)
w * Permutation([2, 1, 3])     # composition (group product)
~w                             # not used; use w.inverse() if needed
w(1)                           # apply to a point: returns w[0]
w[0]                           # 0-indexed array access
len(w)                         # length of the underlying array
list(w)                        # explicit array form
```

### Construction shortcuts

| Call | Meaning |
|---|---|
| `Permutation([3, 1, 2])` | from one-line |
| `Permutation.from_code([2, 0])` or `uncode([2, 0])` | from Lehmer code |
| `Permutation.from_cycles([(1, 3), (2, 4)])` | from cycle notation |
| `Permutation.w0(n)` | longest element of $S_n$ |
| `Permutation.all_permutations(n)` | iterator over $S_n$ |
| `Permutation.longest_element(*descs)` | longest element of a parabolic |
| `Permutation.reflection(root)` | simple/general reflection |

### Useful methods (selection)

`apply`, `swap`, `dominates`, `bruhat_leq`, `min_coset_rep`, `max_coset_rep`,
`coset_decomp`, `parabolic_reduce`, `weight_coset_decomp`,
`min_of_weight_coset`, `max_of_weight_coset`, `code_word`, `right_act`,
`act_root`, `root_swap`, `pad_code`, `mul_dominant`, `strict_mul_dominant`,
`pivots`, `rothe_diagram`, `diagram`, `has_pattern`, `get_cycles`,
`from_partial`, `cycle`.

The full source is [`combinatorics/permutation.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/combinatorics/permutation.py).

---

## Schubert rings

The recommended high-level interface. Elements behave like dicts
`{Permutation: coefficient}` but expose ring operations.

```python
from schubmult import Sx, DSx, Permutation
```

### `Sx` — single (ordinary) Schubert ring

```python
a = Sx([3, 1, 2])              # S_{312}
b = Sx([2, 1, 3])              # S_{213}
c = a * b                      # multiplication via schubmult_py kernel
c                              # {Permutation([4, 1, 2, 3]): 1}
(a + b).expand()               # explicit polynomial in x_1, x_2, ...
list(c.items())                # iterate (Permutation, coeff) pairs
```

`Sx(arg)` accepts a list/tuple (one-line), a `Permutation`, or a polynomial
expression in the underlying variable set.

### `DSx` — double Schubert ring

```python
from schubmult import DSx, Permutation
a = DSx([3, 1, 2])             # double Schubert in y-variables
b = DSx([2, 1, 3])
c = a * b                      # coefficients are polynomials in y
```

For mixed-variable products (`y` *and* `z`), instantiate the ring with two
generating sets — see [`rings/schubert/double_schubert_ring.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/rings/schubert/double_schubert_ring.py).

### Quantum and parabolic-quantum rings

Available via lazy import; pick the right one for your context:

```python
from schubmult import (
    SingleSchubertRing, DoubleSchubertRing,            # ordinary
    # Quantum variants live in:
    #   schubmult.rings.schubert.quantum_schubert_ring
    #   schubmult.rings.schubert.quantum_double_schubert_ring
    #   schubmult.rings.schubert.parabolic_quantum_schubert_ring
    #   schubmult.rings.schubert.parabolic_quantum_double_schubert_ring
)
```

For a custom variable set:

```python
from schubmult import SingleSchubertRing, GeneratingSet
R = SingleSchubertRing(GeneratingSet("u"))      # Schubert ring in u_1, u_2, ...
```

### Common ring-element operations

| Operation | Meaning |
|---|---|
| `a + b`, `a - b`, `a * b`, `n * a` | ring arithmetic |
| `a == b` | structural equality (same dict) |
| `a.expand()` | explicit polynomial in the underlying variables |
| `dict(a)` / `a.items()` | inspect the `{Permutation: coeff}` representation |
| `R.zero`, `R.one` | additive/multiplicative identities of ring `R` |
| `R.new(x)` | coerce a `Permutation`, list, tuple, or expression into `R` |
| `R.from_dict({...})` | build an element directly from a dict |

---

## Low-level multiplication kernels

When you don't want a ring wrapper, work with raw dicts.

### `schubmult_py(perm_dict, v)`

Ordinary (single) Schubert multiplication.
[`mult/single.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/mult/single.py)

```python
from schubmult.mult.single import schubmult_py
from schubmult import Permutation
from sympy import S

u, v = Permutation([3, 1, 2]), Permutation([2, 1, 3])
schubmult_py({u: S.One}, v)
# {Permutation([4, 1, 2, 3]): 1}
```

### `schubmult_double(perm_dict, v, var2=None, var3=None)`

Double Schubert multiplication; coefficients are polynomials in `var2`
(and optionally `var3` for mixed-variable). [`mult/double.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/mult/double.py)

Convenience helpers:

- `schubmult_double_pair(perm1, perm2, var2, var3)` — multiply two Schubert
  basis elements directly, returning the coefficient dict.
- `schubmult_double_pair_generic(perm1, perm2)` — same with the generic
  variable set.
- `mult_poly_double(coeff_dict, poly)` — multiply a Schubert expansion by an
  ordinary polynomial.

### `schubmult_q(perm_dict, v)` and `schubmult_q_fast(...)`

Quantum Schubert multiplication; coefficients live in
$\mathbb{Z}[q_1, q_2, \ldots]$. [`mult/quantum.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/mult/quantum.py)

### `schubmult_q_double(...)`

Quantum double; conjectural in some cases. [`mult/quantum_double.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/mult/quantum_double.py)

### `schub_coprod_py(perm, indices)`

Coproduct of an ordinary Schubert polynomial along the split `indices`.

### Dict utilities

```python
from schubmult.utils.perm_utils import add_perm_dict
add_perm_dict(d1, d2)          # merge two coefficient dicts
```

---

## RC-graphs and root tableaux

Combinatorial models of reduced expressions used internally and exposed for
research.

```python
from schubmult import RCGraph, RootTableau
g = RCGraph.from_perm(Permutation([3, 1, 2]))
g.normalize()                  # canonical form
g.weight                       # weight vector (tuple of ints)
g.raising_operator(1)          # crystal raising op; returns None if undefined
g.lowering_operator(1)
g.to_highest_weight()          # (highest-weight element, sequence of ops)
```

Crystal-operator results are `None` when undefined — always check before
using.

`RootTableau` implements dual Knuth equivalence with JDT slides and tracks
the Edelman–Greene invariant; see
[`combinatorics/root_tableau.py`](https://github.com/matthematics/schubmult/blob/main/src/schubmult/combinatorics/root_tableau.py).

---

## Symbolic layer

`schubmult.symbolic` is a thin facade over SymPy + SymEngine. **Always import
through it**, not directly from sympy/symengine, so the package can swap
backends transparently.

```python
from schubmult.symbolic import sympify, Add, Mul, S
from schubmult.symbolic.poly.variables import GeneratingSet

X = GeneratingSet("x")         # x_1, x_2, ... on demand
X[1] + X[2]                    # symbolic expression
sympify("y_1 + y_2")
```

`GeneratingSet(label)` returns an indexed family; the existing rings expose
their generating set via `R.genset` (and `R.coeff_genset` for double rings).

---

## CLI scripts

Each is registered as a console entry point and also importable as a module
under `schubmult._scripts`.

| Command | Module | Purpose |
|---|---|---|
| `schubmult_py` | `schubmult._scripts.schubmult_py` | Ordinary Schubert product |
| `schubmult_double` | `schubmult._scripts.schubmult_double` | Double Schubert product |
| `schubmult_q` | `schubmult._scripts.schubmult_q` | Quantum Schubert product |
| `schubmult_q_double` | `schubmult._scripts.schubmult_q_double` | Quantum double Schubert product |

Common flags (see `schubmult/utils/argparse.py`):

| Flag | Meaning | Available on |
|---|---|---|
| `--code` | Treat tokens as Lehmer codes instead of permutations | all |
| `-np`, `--no-print` | Suppress stdout | all |
| `--mult EXPR` | Multiply by a polynomial parsed via SymPy | all |
| `--display-positive` | MILP-based root-positive display (slow) | double, q_double |
| `--mixed-var` | Use `y` and `z` instead of just `y` | double, q_double |
| `--down` | Compute the "down" version (divided-difference style) | double |
| `--coprod` | Coproduct mode: one perm + split-index list | py, double |
| `--parabolic g1 g2 …` | Parabolic generators (positive integers) | q, q_double |
| `--basic-pieri` | Slow reference Pieri implementation | q, q_double |
| `--display-mode {basic,pretty,latex,raw}` | Output formatting | quantum scripts |

Examples (terminal):

```bash
schubmult_py 3 1 2 - 2 1 3
schubmult_double 3 1 2 - 2 1 3 --mixed-var --display-positive
schubmult_q 3 1 2 - 2 1 3 --parabolic 1
schubmult_py --code 2 0 - 1 0
```

Programmatic use:

```python
from schubmult._scripts import schubmult_py
schubmult_py.main(["schubmult_py", "3", "1", "2", "-", "2", "1", "3"])
# (prints to stdout; capture with contextlib.redirect_stdout if needed)
```

---

## Worked example: a one-page tour

```python
from schubmult import Sx, DSx, Permutation, uncode

# 1. Build a Schubert polynomial.
a = Sx([3, 1, 2])
print(a.expand())                       # x_1 * x_2

# 2. Multiply two of them.
b = Sx([2, 1, 3])
c = a * b                               # {Permutation([4,1,2,3]): 1}
for w, coeff in c.items():
    print(coeff, w.trimcode)

# 3. Same product with double Schuberts (coefficients in y).
A = DSx([3, 1, 2])
B = DSx([2, 1, 3])
C = A * B                               # coefficients are polynomials in y

# 4. Permutation gymnastics.
w = uncode([2, 0, 1])                   # Lehmer code -> Permutation
print(w, w.code, w.descents())
print(w * Permutation([2, 1, 3]))       # group multiplication

# 5. Hit a kernel directly for raw dict output.
from schubmult.mult.single import schubmult_py
from sympy import S
print(schubmult_py({Permutation([3,1,2]): S.One}, Permutation([2,1,3])))
```
