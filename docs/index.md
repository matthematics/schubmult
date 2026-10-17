# schubmult

**Fast Schubert calculus in Python** — Littlewood-Richardson coefficients for Schubert, double,
quantum, quantum double and Grothendieck polynomials; the polynomial ring in its combinatorial bases;
and the objects that index all of it (permutations, pipe dreams, bumpless pipe dreams, tableaux and
their crystals). Products run on dedicated kernels built on SymEngine, so results are fast and drop
straight into SymPy or SageMath.

[:material-download: PyPI](https://pypi.org/project/schubmult/){ .md-button }
[:material-github: Source](https://github.com/matthematics/schubmult){ .md-button }
[:material-book-open-variant: Module reference](modules/README.md){ .md-button .md-button--primary }

---

## Install

```bash
pip install schubmult             # latest release (Python 3.10+)
pip install --pre schubmult       # latest pre-release
sage -pip install --pre schubmult # inside SageMath, for schubmult.sage
```

SageMath is optional: everything except `schubmult.sage` works in plain Python.

## Five-minute tour

Rings are indexed by permutations in one-line notation; elements are `{Permutation: coefficient}`
dictionaries with `*`, `+`, `.expand()`, `.coproduct()` and change of basis. Variables are
1-indexed on the Python side (`x_1, x_2, ...`).

```python
from schubmult import Sx, DSx, QSx, QDSx, Gx, DGx, Permutation
from schubmult.abc import x, y, z

Sx([3, 1, 2]) * Sx([1, 3, 2])                 # S_312 · S_132
# Sx((3, 2, 1)) + Sx((4, 1, 2, 3))

Sx.from_expr(x[1]**2 + x[1] * x[2])           # any polynomial back into the Schubert basis
# Sx((2, 3, 1)) + Sx((3, 1, 2))

DSx([1, 3, 2]) * DSx([2, 1, 3], z)            # double Schubert, mixed alphabets S_w(x; y) · S_v(x; z)
# (y_1 - z_1)*DSx((1, 3, 2), y) + DSx((2, 3, 1), y) + DSx((3, 1, 2), y)

QSx([2, 1, 3]) * QSx([2, 1, 3])               # quantum Monk: S_{s_1}^2 = S_312 + q_1
# q_1 + QSx((3, 1, 2))

(DGx([2, 1, 3]) * DGx([2, 1, 3])).simplify()  # equivariant K-theory: rational in y and β
# (y_1 - y_2)/(1 + β*y_2)*DGx((2, 1), y) + (1 + β*y_1)/(1 + β*y_2)*DGx((3, 1, 2), y)

w = Permutation([3, 1, 4, 2])
w.code, w.inv, w.descents()                   # ([2, 0, 1], 3, {0, 2})
```

!!! tip "Hidden zeros"
    Kernels leave coefficients as products of factors, and in big products some of them cancel to
    zero without looking like it. `elem.strip_zeros(exact=True)` removes them at a fraction of the
    cost of expanding: `QDSx([2, 4, 1, 3]) * QDSx([3, 1, 4, 2])` has 11 terms, 10 of them nonzero.

The rings and what they compute:

| ring | polynomials | coefficients |
|---|---|---|
| `Sx` | Schubert $\mathfrak S_w(x)$ | integers |
| `DSx` | double Schubert $\mathfrak S_w(x; y)$ | polynomials in $y$ (and $z$ for mixed products) |
| `QSx`, `QPSx(*blocks)` | quantum Schubert, parabolic quantum | polynomials in $q$ |
| `QDSx`, `QPDSx(*blocks)` | quantum double, parabolic quantum double | polynomials in $q, y$ |
| `Gx` | $\beta$-Grothendieck $\mathfrak G^\beta_w(x)$ | polynomials in $\beta$ |
| `DGx` | double Grothendieck $\mathfrak G^\beta_w(x; y)$ | rational functions in $\beta, y$ |

See [`schubmult.rings.schubert`](modules/schubmult/rings/schubert/index.md) for the ring classes
and [`schubmult.mult`](modules/schubmult/mult/index.md) for the multiplication kernels underneath
([`single`](modules/schubmult/mult/single.md), [`double`](modules/schubmult/mult/double.md),
[`quantum`](modules/schubmult/mult/quantum.md), [`quantum_double`](modules/schubmult/mult/quantum_double.md),
[`groth`](modules/schubmult/mult/groth.md), [`groth_double`](modules/schubmult/mult/groth_double.md)).

## Where to look

| I want to... | start here |
|---|---|
| multiply Schubert polynomials of any flavor | [`rings.schubert`](modules/schubmult/rings/schubert/index.md) — `Sx`, `DSx`, `QSx`, `QDSx`, `Gx`, `DGx` |
| work with permutations, Lehmer codes, Bruhat order | [`combinatorics.permutation`](modules/schubmult/combinatorics/permutation.md) |
| enumerate pipe dreams / RC graphs and apply crystal operators | [`combinatorics.rc_graph`](modules/schubmult/combinatorics/rc_graph.md), [`crystal_graph`](modules/schubmult/combinatorics/crystal_graph.md) |
| bumpless pipe dreams and the Gao-Huang bijection | [`combinatorics.bpd`](modules/schubmult/combinatorics/bpd.md), [`mbpd`](modules/schubmult/combinatorics/mbpd.md) |
| tableaux: plactic, nil-plactic, root, increasing, set-valued | [`plactic`](modules/schubmult/combinatorics/plactic.md), [`nilplactic`](modules/schubmult/combinatorics/nilplactic.md), [`root_tableau`](modules/schubmult/combinatorics/root_tableau.md), [`increasing_tableau`](modules/schubmult/combinatorics/increasing_tableau.md), [`set_valued_tableau`](modules/schubmult/combinatorics/set_valued_tableau.md) |
| expand polynomials in the key, slide, forest, glide, Lascoux, grove bases | [`rings.polynomial_algebra`](modules/schubmult/rings/polynomial_algebra/index.md) |
| the dual picture: the free algebra and its bases | [`rings.free_algebra`](modules/schubmult/rings/free_algebra/index.md) |
| write double coefficients as manifestly positive expressions | [`mult.positivity`](modules/schubmult/mult/positivity.md), or `--display-positive` on the command line |
| use all of this inside SageMath | [`schubmult.sage`](modules/schubmult/sage/index.md) — see below |
| the command-line tools | [`schubmult._scripts`](modules/schubmult/_scripts/index.md) |

## Combinatorics

The [`schubmult.combinatorics`](modules/schubmult/combinatorics/index.md) package provides the
objects that index Schubert calculus, with conversions between them:

```python
from schubmult import RCGraph, BPD, Permutation

w = Permutation([1, 4, 2, 3])
rcs = RCGraph.all_rc_graphs(w, 3)        # the 3 RC graphs of w in 3 rows

rc = next(rc for rc in rcs if rc.lowering_operator(1) is not None)
rc.lowering_operator(1)                  # Kashiwara operator f_1 (None when undefined)

bpd = BPD.from_rc_graph(rc)              # Gao-Huang bijection to bumpless pipe dreams
assert bpd.to_rc_graph() == rc
```

Crystal operators return `None` where they are undefined rather than raising; highest weights come
from `.to_highest_weight()`.

## SageMath

[`schubmult.sage`](modules/schubmult/sage/index.md) provides Sage parents (built on
`CombinatorialFreeModule`, with coercion, `TestSuite`-clean) whose arithmetic is delegated to the
kernels. Variables are 0-indexed there, as in Sage (`x0`, `y_0`, `q_0`).

```python
sage: from schubmult.sage import *
sage: X = DoubleSchubertPolynomialRing(QQ)
sage: X([3, 1, 2]) * X([2, 1])
(y_2-y_0)*X_y[3, 1, 2] + X_y[4, 1, 2, 3]

sage: QD = QuantumDoubleSchubertPolynomialRing(QQ, parabolic=(2, 3))   # QH_T^*(Gr(2, 5))
sage: (QD([2, 4, 1, 3]) * QD([2, 5, 1, 3, 4])).project(5)               # Buch's [X^(2,1)] * [X^(3,1)]
(q_0*y_4*y_1-q_0*y_4*y_0-q_0*y_2*y_1+q_0*y_2*y_0)*Xq_y[1] + (q_0*y_4-q_0*y_0)*Xq_y[1, 3, 2] + q_0*Xq_y[1, 4, 2, 3] + q_0*Xq_y[2, 3, 1] + (y_4^2*y_1-y_4^2*y_0-y_4*y_2*y_1+y_4*y_2*y_0-y_4*y_1*y_0+y_4*y_0^2+y_2*y_1*y_0-y_2*y_0^2)*Xq_y[2, 5, 1, 3, 4] + (y_4^2-2*y_4*y_0+y_0^2)*Xq_y[3, 5, 1, 2, 4] + (y_4-y_0)*Xq_y[4, 5, 1, 2, 3]

sage: GD = DoubleGrothendieckPolynomialRing(QQ)                         # K_T of the flag variety
sage: GD([2, 1]) * GD([2, 1])
-((y_1-y_0)/(beta*y_1+1))*G_y[2, 1] + ((beta*y_0+1)/(beta*y_1+1))*G_y[3, 1, 2]

sage: A = PolynomialAlgebra(QQ); S = A.schubert(); k = A.key()         # QQ[x0, x1, ...] with its bases
sage: k(S[3, 1, 2])
k[2]
sage: A.fundamental_slide()(KeyPolynomials(QQ)([1, 0, 2]))
F[1, 0, 2] + F[2, 0, 1]
```

- [`DoubleSchubertPolynomialRing`](modules/schubmult/sage/double_schubert.md),
  [`QuantumSchubertPolynomialRing` / `QuantumDoubleSchubertPolynomialRing`](modules/schubmult/sage/quantum_schubert.md)
  (with `parabolic=` block sizes), [`GrothendieckPolynomialRing` / `DoubleGrothendieckPolynomialRing`](modules/schubmult/sage/grothendieck.md).
- [`PolynomialAlgebra(R)`](modules/schubmult/sage/polynomial_algebra.md) and `PolynomialAlgebra(R, n)`:
  the polynomial ring with its bases as mutually coercing realizations — Schubert and Grothendieck
  indexed by permutations; monomial, key, slide, forest, glide, Lascoux, grove by weak compositions;
  and in `n` variables the elementary symmetric basis.
- Sage's own `SchubertPolynomialRing`, `KeyPolynomials` and `SymmetricFunctions` coerce or convert in
  and out (`from_symmetric_function`, `to_symmetric_function`).

!!! warning "Conventions"
    - **β.** The Grothendieck *rings* keep `beta` as a variable. The K-theoretic bases of
      `PolynomialAlgebra` (`grothendieck`, `glide`, `lascoux`, `grove`) are at **β = −1**, the
      classical convention: `G[1, 3, 2].expand()` is `x0 + x1 - x0*x1`.
    - **Double Grothendieck** polynomials use $x \oplus y = x + y + \beta x y$, so at β = 0 they are
      double Schubert polynomials in $x$ and $-y$.
    - **Stability.** The rings are stable, so a product carries every class of the infinite flag
      variety; `project(n)` drops the classes that vanish in the flag variety of $\mathbb C^n$.

## Command line

Each ring has a script; permutations are space-separated and factors separated by `-`:

```bash
schubmult_py 3 1 2 - 2 1 3                         # Schubert
schubmult_py --code 2 0 - 1 0                      # the same product via Lehmer codes
schubmult_double 1 3 2 - 1 3 2 --display-positive  # double, coefficients written Graham-positively
schubmult_q 2 1 3 - 2 1 3                          # quantum
schubmult_q_double 2 1 3 - 2 1 3 --parabolic 1     # parabolic quantum double
grothmult_py 2 1 3 - 2 1 3                         # Grothendieck
grothmult_double 2 1 3 - 2 1 3                     # double Grothendieck
```

Run any script with `--help` for the full option list; see [`schubmult._scripts`](modules/schubmult/_scripts/index.md).

## Contributing

Issues and pull requests are welcome at [github.com/matthematics/schubmult](https://github.com/matthematics/schubmult).
Feature branches target `develop`; `main` carries releases. The
[CHANGELOG](https://github.com/matthematics/schubmult/blob/main/CHANGELOG.md) lists what changed in
each version. schubmult is GPL-3.0 licensed.
