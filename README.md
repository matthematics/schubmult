# schubmult

**Fast Schubert calculus in Python: Littlewood-Richardson coefficients, Schubert/Grothendieck rings, and the combinatorics behind them.**

[![PyPI](https://img.shields.io/pypi/v/schubmult)](https://pypi.org/project/schubmult/)
[![Docs](https://img.shields.io/badge/docs-matthematics.github.io%2Fschubmult-blue)](https://matthematics.github.io/schubmult/)
[![License: GPL v3](https://img.shields.io/badge/license-GPLv3-blue.svg)](LICENSE)

`schubmult` multiplies single, double, quantum, and quantum double Schubert polynomials (with parabolic variants), expands products in the Schubert basis, and gives you the combinatorial objects that index them: permutations, RC graphs (pipe dreams), bumpless pipe dreams, tableaux, and their crystal structures. It is built on [SymEngine](https://github.com/symengine/symengine.py) for speed and [SymPy](https://www.sympy.org/) for display, so results drop straight into SymPy or Sage.

**Documentation:** <https://matthematics.github.io/schubmult/>

---

## Installation

```bash
pip install schubmult
```

Development version:

```bash
pip install git+https://github.com/matthematics/schubmult.git
```

Requires Python 3.10+. To build the documentation locally, install the `docs` extra: `pip install -e ".[docs]"`.

## Quick start

### Python

```python
from schubmult import Sx, DSx, QSx, QDSx, Gx, DGx, Permutation
from schubmult.abc import x, y, z

# Ordinary Schubert polynomials: S_312 * S_132 = S_321 + S_4123
Sx([3, 1, 2]) * Sx([1, 3, 2])
# 𝔖_(3, 2, 1)(x) + 𝔖_(4, 1, 2, 3)(x)

# Expand to a polynomial, or go the other way
Sx([3, 1, 2]).expand()                       # x_1**2
Sx.from_expr(x[1]**2 + x[1] * x[2])          # 𝔖_(2, 3, 1)(x) + 𝔖_(3, 1, 2)(x)

# Double Schubert polynomials S_w(x; y); coefficients are polynomials in y
DSx([1, 3, 2]) * DSx([1, 3, 2])
# (-y_2 + y_3)*𝔖_(1, 3, 2)(x; y) + 𝔖_(1, 4, 2, 3)(x; y) + 𝔖_(2, 3, 1)(x; y)

# Mixed variables: S_w(x; y) * S_v(x; z)
DSx([1, 3, 2]) * DSx([2, 1, 3], z)
# (y_1 - z_1)*𝔖_(1, 3, 2)(x; y) + 𝔖_(2, 3, 1)(x; y) + 𝔖_(3, 1, 2)(x; y)

# Quantum Schubert polynomials
QSx([2, 1, 3]) * QSx([2, 1, 3])              # q_1 + 𝕼𝔖_312(x)

# Quantum double Schubert polynomials
QDSx([2, 1, 3]) * QDSx([2, 1, 3])
# q_1 + (-y_1 + y_2)*𝕼𝔖_21(x; y) + 𝕼𝔖_312(x; y)

# Grothendieck polynomials (K-theory)
Gx([2, 1, 3]) * Gx([2, 1, 3])                # 𝔊_(3, 1, 2)(x)

# Double Grothendieck polynomials; .simplify() puts the rational coefficients in y, β in normal form
(DGx([1, 3, 2]) * DGx([2, 1, 3])).simplify()
# β*𝔊_(3, 2, 1)(x; y) + 𝔊_(2, 3, 1)(x; y) + 𝔊_(3, 1, 2)(x; y)
(DGx([2, 1, 3]) * DGx([2, 1, 3])).simplify()
# (y_1 - y_2)*𝔊_(2, 1)(x; y)   (y_1*β + 1)*𝔊_(3, 1, 2)(x; y)
# -------------------------- + -----------------------------
#         y_2*β + 1                      y_2*β + 1

# Permutations use one-line notation; Lehmer codes are available too
w = Permutation([3, 1, 4, 2])
w.code, w.inv, w.descents()                  # ([2, 0, 1], 3, {0, 2})
```

Ring objects: `Sx` (single), `DSx` (double), `QSx`/`QDSx` (quantum, quantum double), `QPSx`/`QPDSx` (parabolic quantum), `Gx`/`DGx` (single/double Grothendieck). Elements are dictionaries `{Permutation: coefficient}` with `*`, `+`, `.expand()`, `.coproduct()`, and change of basis between them. The underlying multiplication kernels live in `schubmult.mult` (`single`, `double`, `quantum`, `quantum_double`, `groth`, `groth_double`).

### Command line

Each ring has a CLI. Permutations are space-separated, factors are separated by `-`:

```bash
schubmult_py 3 1 2 - 2 1 3                         # 1  (4, 1, 2, 3)
schubmult_py --code 2 0 - 1 0                      # same product via Lehmer codes
schubmult_double 1 3 2 - 1 3 2 --display-positive  # double Schubert, coefficients displayed positively
schubmult_q 2 1 3 - 2 1 3                          # quantum
schubmult_q_double 2 1 3 - 2 1 3 --parabolic 1     # parabolic quantum double
grothmult_py 2 1 3 - 2 1 3                         # Grothendieck
grothmult_double 2 1 3 - 2 1 3                     # double Grothendieck (coefficients in y and β)
```

`--display-positive` writes double and quantum double coefficients as manifestly positive expressions in the differences `y_i - z_j` (Graham positivity), using integer programming to find a positive representative. Run any script with `--help` for the full option list.

## Combinatorics

The `schubmult.combinatorics` package provides the objects that index Schubert calculus, with conversions between them:

- **`Permutation`** -- one-line notation, Lehmer codes (`uncode`), reduced words, Bruhat order, descents, dominant/Grassmannian tests.
- **`RCGraph`** -- reduced compatible sequences / pipe dreams, with enumeration (`RCGraph.all_rc_graphs(w, n)`), Kashiwara crystal operators (`raising_operator(i)`, `lowering_operator(i)`), Edelman-Greene insertion, and the crystal products used in the transition formulas.
- **`BPD`** -- bumpless pipe dreams with droop moves, the Gao-Huang bijection to RC graphs (`BPD.from_rc_graph`, `.to_rc_graph()`), and marked/unreduced variants for Grothendieck polynomials.
- **`WCGraph`**, **`PipeDream`**, **`HPD`** -- K-theoretic and hybrid pipe dream models.
- **Tableaux** -- `Plactic` (semistandard, jeu de taquin), `NilPlactic` (Edelman-Greene), `RootTableau`, `IncreasingTableau`, `SetValuedTableau`, `HeckePlactic`.
- **Forests** -- indexed forests and the Thompson monoid factorization behind the forest and grove bases.

```python
from schubmult import RCGraph, BPD, Permutation

w = Permutation([1, 4, 2, 3])
rcs = RCGraph.all_rc_graphs(w, 3)        # all RC graphs of w in 3 rows

rc = next(rc for rc in rcs if rc.lowering_operator(1) is not None)
rc.lowering_operator(1)                  # crystal operator f_1 (None when undefined)

bpd = BPD.from_rc_graph(rc)              # Gao-Huang bijection
assert bpd.to_rc_graph() == rc
```

## Rings and algebras

Beyond the Schubert rings, `schubmult.rings` includes:

- **`PolynomialAlgebra`** with interchangeable bases: monomials, Schubert, key (Demazure), Lascoux, fundamental/monomial slide, glide, forest, grove, Grothendieck, and elementary symmetric bases.
- **`FreeAlgebra`** -- the graded dual of the polynomial algebra (a word `(a_1, ..., a_n)` is dual to `x_1^a_1 ... x_n^a_n`), with a dual basis for each of the above; `ASx` is the dual Schubert basis.
- **Combinatorial rings** -- `RCGraphRing`, `WCGraphRing`, `CrystalGraphRing`, and related algebras whose basis elements are the combinatorial objects themselves.
- **Tensor and direct products**, `NSym`, `QSym`, and the nil-Hecke algebra.

```python
from schubmult import ASx, FA

ASx([2, 1, 3]) * ASx([1, 3, 2])   # dual Schubert basis of the free algebra
FA(1, 0) * FA(2)                  # word basis: concatenation, [102]
```

## Symbolic layer

`schubmult.symbolic` wraps SymEngine and SymPy behind one interface (`sympify`, `expand`, `Add`, `Mul`, `S`), provides indexed variable families (`GeneratingSet("x")` gives `x_1, x_2, ...`), and unevaluated elementary/complete symmetric polynomial atoms (`E`, `e`, `H`, `h`) so Schubert polynomials can be manipulated in the SEM basis. `schubmult.abc` exposes ready-made `x`, `y`, `z`, `q`, `beta` in the spirit of `sympy.abc`.

## Documentation

Full API documentation, generated from the source docstrings, is at **<https://matthematics.github.io/schubmult/>**. To build it locally:

```bash
pip install -e ".[docs]"
python docs/generate_docs.py
mkdocs serve
```

## Development

```bash
git clone https://github.com/matthematics/schubmult.git
cd schubmult
pip install -e .
pytest
```

Tests live in `tests/`; the script tests compare CLI output against stored JSON cases in `tests/scripts/data`. `ruff check src/schubmult` should pass.

## License

GPL-3.0. See [LICENSE](LICENSE).
