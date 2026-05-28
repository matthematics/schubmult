"""Compute double forest polynomials P_F(x; t) via the vine subword model.

Reference
---------
N. Bergeron, L. Gagnon, P. Nadeau, H. Spink, V. Tewari,
"Equivariant quasisymmetry and noncrossing partitions" (arXiv:2504.15234),
Theorem 5.1 and Section 5.1.

Given a code c = (c_1, c_2, ...) for an indexed forest F (so that c(F) = c),
and two sequences of symbols x = (x_1, x_2, ...) and t = (t_1, t_2, ...),
this returns the double forest polynomial

    P_F(x; t) = sum_{pi in R(omega_tilde_[n]; F)} wt(pi),

where omega_tilde_[n] is the long vine word and the sum is over its
subwords whose value-sequence is a Sylvester word for F.
"""

from __future__ import annotations

from itertools import combinations
from typing import Dict, Iterable, Tuple

from sympy import IndexedBase, expand, factor


# --------------------------------------------------------------------------- #
# Forest representation                                                       #
# --------------------------------------------------------------------------- #
# A binary tree is either
#   ("leaf", n)        -- a leaf with absolute label n in N, or
#   ("node", L, R)     -- an internal node with binary-tree children L, R.
# A forest is a tuple of binary trees, ordered by leaf positions.


def _rightmost_leaf(t):
    while t[0] == "node":
        t = t[2]
    return t[1]


def _max_leaf(t):
    if t[0] == "leaf":
        return t[1]
    return max(_max_leaf(t[1]), _max_leaf(t[2]))


def _mul_elementary(forest, i):
    """Right-multiply ``forest`` by the elementary forest ``i_`` in the
    Thompson monoid: leaves > i shift up by 1, and leaf i is wrapped by a
    new internal node whose right child is a brand-new leaf i+1."""

    def shift(t):
        if t[0] == "leaf":
            return ("leaf", t[1] + 1) if t[1] > i else t
        return ("node", shift(t[1]), shift(t[2]))

    forest = [shift(t) for t in forest]
    found = [False]

    def wrap(t):
        if found[0]:
            return t
        if t[0] == "leaf":
            if t[1] == i:
                found[0] = True
                return ("node", ("leaf", i), ("leaf", i + 1))
            return t
        nl = wrap(t[1])
        if found[0]:
            return ("node", nl, t[2])
        return ("node", t[1], wrap(t[2]))

    forest = [wrap(t) for t in forest]
    if not found[0]:
        existing = max((_max_leaf(t) for t in forest), default=0)
        for k in range(existing + 1, i):
            forest.append(("leaf", k))
        forest.append(("node", ("leaf", i), ("leaf", i + 1)))
    return forest


def forest_from_code(code):
    """Build the indexed forest F with c(F) = ``code`` via the Thompson
    monoid factorization F = 1^{c_1} . 2^{c_2} . 3^{c_3} . ..."""
    F = []
    for i, c in enumerate(code, start=1):
        for _ in range(c):
            F = _mul_elementary(F, i)
    return tuple(F)


def sylvester_word(forest):
    """Return one Sylvester word of ``forest`` via pre-order traversal
    (root, left subtree, right subtree) over each non-trivial tree."""
    word = []

    def pre(t):
        if t[0] == "leaf":
            return
        word.append(_rightmost_leaf(t[1]))  # canonical label
        pre(t[1])
        pre(t[2])

    for t in forest:
        pre(t)
    return word


# --------------------------------------------------------------------------- #
# Canonical form for testing whether two words are Sylvester-equivalent       #
# --------------------------------------------------------------------------- #
# Every injective word w is a Sylvester word of a unique indexed forest,
# obtained as follows (Section 3.1 of the paper):
#   - Partition supp(w) into maximal contiguous integer intervals.
#   - For each interval I, take the subword of w using only letters in I.
#   - The BST built by inserting that subword left-to-right is the labeled
#     internal-node skeleton of one tree of the forest.
# Two words are Sylvester-equivalent iff their intervals and per-interval
# BSTs coincide.


def _bst_insert(t, a):
    if t is None:
        return (a, None, None)
    r, L, R = t
    if a < r:
        return (r, _bst_insert(L, a), R)
    return (r, L, _bst_insert(R, a))


def canonical_forest_from_word(word):
    if not word:
        return ()
    sup = sorted(set(word))
    intervals = []
    cur = [sup[0]]
    for x in sup[1:]:
        if x == cur[-1] + 1:
            cur.append(x)
        else:
            intervals.append((cur[0], cur[-1]))
            cur = [x]
    intervals.append((cur[0], cur[-1]))
    out = []
    for lo, hi in intervals:
        iset = set(range(lo, hi + 1))
        sub = [a for a in word if a in iset]
        bst = None
        for a in sub:
            bst = _bst_insert(bst, a)
        out.append((lo, hi, bst))
    return tuple(out)


# --------------------------------------------------------------------------- #
# The vine long word and weights                                              #
# --------------------------------------------------------------------------- #


def long_word(n):
    """Build omega_tilde_[n] as a list of letter records.

    omega^(k)_[n] = (n, n-1, ..., k+1, k, k+1_bar, ..., n-1_bar, n_bar)
    Concatenate for k = 1, ..., n.
    """
    letters = []
    for k in range(1, n + 1):
        for v in range(n, k - 1, -1):
            letters.append({"value": v, "barred": False, "syllable": k})
        for v in range(k + 1, n + 1):
            letters.append({"value": v, "barred": True, "syllable": k})
    return letters


def letter_weight(letter, x_gen, t_gen):
    """wt(j^(i))      = x_i - t_j   (unbarred)
    wt(j_bar^(i)) = t_j - t_i   (barred)"""
    i = letter["syllable"]
    j = letter["value"]
    if letter["barred"]:
        return t_gen(j) - t_gen(i)
    return x_gen(i) - t_gen(j)


# --------------------------------------------------------------------------- #
# Main entry point                                                            #
# --------------------------------------------------------------------------- #


def double_forest_polynomial(code, x_gen, t_gen, n=None):
    """Compute P_F(x; t) for the indexed forest F with c(F) = ``code``.

    Parameters
    ----------
    code : sequence[int]
        The code (c_1, c_2, ...) of F. Trailing zeros are ignored.
    x_gen : callable[int -> Expr]
        Maps i -> x_i  (the non-equivariant variables).
    t_gen : callable[int -> Expr]
        Maps j -> t_j  (the equivariant variables).
    n : int, optional
        Use the long word omega_tilde_[n]. If ``None``, the smallest n with
        F supported in {1, ..., n+1} is used.

    Returns
    -------
    sympy.Expr
        The expanded polynomial P_F(x; t).
    """
    code = list(code)
    while code and code[-1] == 0:
        code.pop()
    F = forest_from_code(code)
    target = canonical_forest_from_word(sylvester_word(F))
    sz = sum(code)  # number of internal nodes = |F|

    max_leaf = 0
    for t in F:
        max_leaf = max(max_leaf, _max_leaf(t))
    inferred_n = max(max_leaf - 1, sz, 1)
    if n is None:
        n = inferred_n
    elif n < inferred_n:
        raise ValueError(f"n={n} too small; need n >= {inferred_n}")

    word = long_word(n)
    total = 0
    for idx in combinations(range(len(word)), sz):
        values = [word[i]["value"] for i in idx]
        if canonical_forest_from_word(values) != target:
            continue
        w = 1
        for i in idx:
            w = w * letter_weight(word[i], x_gen, t_gen)
        total = total + w
    return expand(total)


def _padded_code(code, length):
    data = list(code)
    if len(data) < length:
        data = data + [0] * (length - len(data))
    return tuple(data)


def _degree_one_code(i, length=None):
    if i < 1:
        raise ValueError("degree-one index i must be >= 1")
    length = max(length or i, i)
    out = [0] * length
    out[i - 1] = 1
    return tuple(out)


def decompose_double_forest_tensor(poly, x_genset, t_genset, length):
    """Decompose a double forest polynomial into t-forest ⊗ x-forest terms.

    Returns a dict mapping ``(t_code, x_code) -> integer/symbolic coefficient``.
    """
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra

    ForestX = PolynomialAlgebra(ForestPolyBasis(x_genset))
    ForestT = PolynomialAlgebra(ForestPolyBasis(t_genset))

    x_expansion = ForestX.from_expr(expand(poly), length=length)
    tensor = {}
    for x_code, t_coeff_expr in x_expansion.items():
        t_expansion = ForestT.from_expr(expand(t_coeff_expr), length=length)
        for t_code, scalar in t_expansion.items():
            key = (tuple(t_code), tuple(x_code))
            tensor[key] = expand(tensor.get(key, 0) + scalar)
    return tensor


class DoubleForestPolynomialBasis:
    """Concrete basis for double forest polynomials.

    A basis element is indexed by a pair ``(t_code, x_code)`` and represents
    ``Forest_t(t_code) ⊗ Forest_x(x_code)``.
    """

    def __init__(self, x_genset, t_genset, length):
        from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra
        from schubmult.rings.tensor_ring import TensorRing

        self.length = length
        self.x_ring = PolynomialAlgebra(ForestPolyBasis(x_genset))
        self.t_ring = PolynomialAlgebra(ForestPolyBasis(t_genset))
        self.ring = TensorRing(self.t_ring, self.x_ring)

    def basis_key(self, t_code, x_code):
        return (_padded_code(t_code, self.length), _padded_code(x_code, self.length))

    def basis_element(self, t_code, x_code, coeff=1):
        key = self.basis_key(t_code, x_code)
        return self.ring.from_dict({key: coeff})

    def from_expr(self, poly):
        tensor_dct = decompose_double_forest_tensor(poly, self.x_ring.genset, self.t_ring.genset, self.length)
        return self.ring.from_dict(tensor_dct)

    def multiply_basis_elements(self, left_key, right_key):
        left = self.basis_element(*left_key)
        right = self.basis_element(*right_key)
        return left * right


class ForestDoubleElement(dict):
    """Formal element in the abstract double-forest basis.

    Stored as {forest_code: coeff}, representing
        sum_F coeff[F] * DF(F),
    where DF(F) is an abstract basis symbol (not expanded into x/t variables).
    """

    def __init__(self, ring, data=None):
        super().__init__()
        self.ring = ring
        if data:
            for k, v in data.items():
                kk = ring.normalize_code(k)
                vv = expand(v)
                if vv != 0:
                    self[kk] = vv

    def __add__(self, other):
        other = self.ring.coerce(other)
        out = dict(self)
        for k, v in other.items():
            out[k] = expand(out.get(k, 0) + v)
            if out[k] == 0:
                del out[k]
        return ForestDoubleElement(self.ring, out)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.__add__(-1 * other)

    def __rmul__(self, scalar):
        if isinstance(scalar, (int, float)) or hasattr(scalar, "free_symbols"):
            return ForestDoubleElement(self.ring, {k: expand(scalar * v) for k, v in self.items()})
        return NotImplemented

    def __mul__(self, other):
        other = self.ring.coerce(other)
        return self.ring.multiply(self, other)

    def expand_polynomial(self):
        """Expand abstract basis expression to a concrete x/t polynomial."""
        return self.ring.expand(self)

    def coefficient(self, code):
        return expand(self.get(self.ring.normalize_code(code), 0))

    def __repr__(self):
        if not self:
            return "0"
        terms = []
        for code, coeff in sorted(self.items()):
            terms.append(f"({coeff})*DF{code}")
        return " + ".join(terms)


class ForestDouble:
    """Abstract double-forest basis algebra.

    Basis symbols are indexed by one forest code F and represent P_F(x; t)
    abstractly. Coefficients live in the polynomial ring of equivariant vars.
    """

    def __init__(self, x_gen, t_gen, x_genset, length, n=None):
        self.x_gen = x_gen
        self.t_gen = t_gen
        self.x_genset = x_genset
        self.length = length
        self.n = n

    def normalize_code(self, code):
        return _padded_code(code, self.length)

    def coerce(self, other):
        if isinstance(other, ForestDoubleElement):
            if other.ring is not self:
                raise TypeError("Cannot combine ForestDoubleElements from different rings")
            return other
        if other == 0:
            return ForestDoubleElement(self, {})
        raise TypeError("Expected ForestDoubleElement")

    def basis(self, code, coeff=1):
        return ForestDoubleElement(self, {self.normalize_code(code): coeff})

    def __call__(self, *code):
        if len(code) == 1 and isinstance(code[0], (tuple, list)):
            return self.basis(code[0])
        return self.basis(code)

    def polynomial_of_basis(self, code):
        return double_forest_polynomial(self.normalize_code(code), self.x_gen, self.t_gen, n=self.n)

    def expand(self, elem: ForestDoubleElement):
        total = 0
        for code, coeff in elem.items():
            total = total + coeff * self.polynomial_of_basis(code)
        return expand(total)

    def from_polynomial(self, poly):
        coeffs = extract_double_forest_coefficients(expand(poly), x_genset=self.x_genset, length=self.length)
        return ForestDoubleElement(self, coeffs)

    def multiply(self, left: ForestDoubleElement, right: ForestDoubleElement):
        # Abstract product is defined by expanding to x/t, multiplying, and
        # re-extracting in the abstract DF basis.
        prod_poly = expand(self.expand(left) * self.expand(right))
        return self.from_polynomial(prod_poly)


def extract_double_forest_coefficients(poly, x_genset, length):
    """Return all coefficients a_F(t) in f = sum_F a_F(t) P_F(x; t).

    Implemented by iterative top-degree peeling:
    1) convert the current residual into ForestPolyBasis in x,
    2) read top-degree coefficients,
    3) subtract those coefficients times the corresponding double forest basis
       polynomials P_F,
    4) repeat until the residual vanishes.

    This avoids treating the double-forest basis as a plain one-shot x-basis
    conversion and follows the triangular subtraction strategy.
    """
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra
    from schubmult.abc import x, y

    ForestX = PolynomialAlgebra(ForestPolyBasis(x_genset))
    residual = {tuple(code): expand(coeff) for code, coeff in ForestX.from_expr(expand(poly), length=length).items() if expand(coeff) != 0}
    out = {}

    # Cache DF(code) expressed in x-forest basis to avoid repeated conversions.
    x_basis_cache = {}

    def df_in_x_basis(code):
        key = tuple(code)
        if key not in x_basis_cache:
            df_poly = double_forest_polynomial(key, lambda i: x[i], lambda j: y[j])
            x_basis_cache[key] = {tuple(k): expand(v) for k, v in ForestX.from_expr(df_poly, length=length).items() if expand(v) != 0}
        return x_basis_cache[key]

    for _ in range(20_000):
        if not residual:
            break

        top_degree = max(sum(code) for code in residual)
        top_codes = [code for code in residual if sum(code) == top_degree]

        progress = False
        for code in top_codes:
            coeff = expand(residual.get(code, 0))
            if coeff == 0:
                continue
            out[code] = expand(out.get(code, 0) + coeff)
            progress = True

            # Subtract coeff * DF(code), but entirely in x-forest-basis dict form.
            for beta, c_beta in df_in_x_basis(code).items():
                residual[beta] = expand(residual.get(beta, 0) - coeff * c_beta)
                if residual[beta] == 0:
                    del residual[beta]

        if not progress:
            raise RuntimeError("Extraction stalled; no progress on top degree terms")

    if residual:
        raise RuntimeError(f"Extraction did not terminate; residual dict has {len(residual)} terms")

    return {code: expand(coeff) for code, coeff in out.items() if expand(coeff) != 0}


def extract_double_forest_coefficient(poly, forest_code, x_genset, length):
    """Return the single coefficient a_F(t) of P_F in f.

    Equivalent to the operator formula a_F(t) = [ev star e_F] f from
    Theorem 10.9 (Coefficient extraction and star-composition).

    Note: this returns a t-polynomial coefficient in the x-forest basis.
    """
    target = _padded_code(forest_code, length)
    coeffs = extract_double_forest_coefficients(poly, x_genset=x_genset, length=length)
    return expand(coeffs.get(target, 0))


def extract_double_forest_tensor_coefficients(poly, x_genset, t_genset, length):
    """Return true double-forest tensor-basis coefficients.

    Expands poly as
        poly = sum_{A,B} c_{A,B} Forest_t(A) \otimes Forest_x(B),
    returning a dict mapping (A, B) -> c_{A,B}.
    """
    return decompose_double_forest_tensor(poly, x_genset=x_genset, t_genset=t_genset, length=length)


def extract_double_forest_tensor_coefficient(poly, t_code, x_code, x_genset, t_genset, length):
    """Return c_{A,B} for one tensor basis pair A=t_code, B=x_code."""
    A = _padded_code(t_code, length)
    B = _padded_code(x_code, length)
    coeffs = extract_double_forest_tensor_coefficients(poly, x_genset=x_genset, t_genset=t_genset, length=length)
    return expand(coeffs.get((A, B), 0))


def monk_style_degree_one_product(code, i, x_gen, t_gen, x_genset, t_genset, n=None, length=None):
    """Compute a Monk-style decomposition for degree-1 double forest multiplication.

    The product
        P_{e_i}(x; t) * P_code(x; t)
    is decomposed as
        diagonal_coeff(t) * P_code(x; t) + sum_{beta != code} c_beta(t) * P_beta(x; t),
    and also expanded in the tensor basis t-forest ⊗ x-forest.
    """
    base_code = tuple(code)
    min_len = max(len(base_code), i)
    if length is None:
        length = min_len
    else:
        length = max(length, min_len)

    left_code = _degree_one_code(i, length=length)
    right_code = _padded_code(base_code, length)

    left_poly = double_forest_polynomial(left_code, x_gen, t_gen, n=n)
    right_poly = double_forest_polynomial(right_code, x_gen, t_gen, n=n)
    product_poly = expand(left_poly * right_poly)

    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra

    ForestX = PolynomialAlgebra(ForestPolyBasis(x_genset))
    x_expansion = {tuple(beta): expand(coeff) for beta, coeff in ForestX.from_expr(product_poly, length=length).items()}

    diagonal_coeff = expand(x_expansion.get(right_code, 0))
    off_diag = {tuple(beta): expand(coeff) for beta, coeff in x_expansion.items() if tuple(beta) != right_code and coeff != 0}
    tensor_expansion = decompose_double_forest_tensor(product_poly, x_genset, t_genset, length=length)
    df_basis = DoubleForestPolynomialBasis(x_genset=x_genset, t_genset=t_genset, length=length)
    double_forest_element = df_basis.ring.from_dict(tensor_expansion)

    return {
        "left_code": left_code,
        "right_code": right_code,
        "left_poly": left_poly,
        "right_poly": right_poly,
        "product_poly": product_poly,
        "diagonal_coeff": diagonal_coeff,
        "off_diagonal": off_diag,
        "double_forest_basis": df_basis,
        "double_forest_element": double_forest_element,
        "tensor_expansion": tensor_expansion,
    }


# --------------------------------------------------------------------------- #
# Self-tests against examples from the paper                                  #
# --------------------------------------------------------------------------- #


def _selftest():
    x = IndexedBase("x")
    t = IndexedBase("t")
    xg = lambda i: x[i]
    tg = lambda j: t[j]

    # Example 4.4: P_2 = x_1 + x_2 - t_1 - t_2
    p2 = double_forest_polynomial([0, 1], xg, tg)
    expected = x[1] + x[2] - t[1] - t[2]
    assert expand(p2 - expected) == 0, f"P_2 mismatch: {p2}"
    print(f"P_2                  = {p2}")

    # Figure 1: F has code (2, 0, 1)
    # P_F = (x_1 - t_1)(x_1 - t_2)(x_2 + x_3 - t_1 - t_2)
    pF = double_forest_polynomial([2, 0, 1], xg, tg)
    expected = expand((x[1] - t[1]) * (x[1] - t[2]) * (x[2] + x[3] - t[1] - t[2]))
    assert expand(pF - expected) == 0, f"P_(2,0,1) mismatch: {pF}"
    print(f"P_(code=(2,0,1))     = {factor(pF)}")

    # Example 4.11: F has code (0, 2, 1) (a zigzag forest)
    # F_(0,2,1) = (x_1-t_4)(x_1-t_1)(x_2+x_3-t_1-t_2)
    #          + (x_1+x_2-t_1-t_4)(x_2-t_2)(x_3-t_2)
    p = double_forest_polynomial([0, 2, 1], xg, tg)
    expected = expand(
        (x[1] - t[4]) * (x[1] - t[1]) * (x[2] + x[3] - t[1] - t[2])
        + (x[1] + x[2] - t[1] - t[4]) * (x[2] - t[2]) * (x[3] - t[2])
    )
    assert expand(p - expected) == 0, f"P_(0,2,1) mismatch: {p}"
    print(f"P_(code=(0,2,1))     = {p}")

    # Spot-check Table 1: P_3 = (x_1-t_3)+(t_3-t_1)+(x_2-t_3)+(t_3-t_2)+(x_3-t_3)
    p3 = double_forest_polynomial([0, 0, 1], xg, tg)
    expected = (
        (x[1] - t[3]) + (t[3] - t[1]) + (x[2] - t[3]) + (t[3] - t[2]) + (x[3] - t[3])
    )
    assert expand(p3 - expected) == 0
    print(f"P_3                  = {p3}")

    # Empty forest: P_empty = 1
    pe = double_forest_polynomial([], xg, tg)
    assert pe == 1
    print(f"P_empty              = {pe}")

    print("\nAll self-tests passed.")


if __name__ == "__main__":
    from sympy import pretty_print
    from symengine import symbols
    from schubmult.abc import x, y
    from schubmult.rings.polynomial_algebra import ForestPolyBasis, PolynomialAlgebra, Forest
    t = [symbols(f"t_{i}") for i in range(1, 10)]
    #from schubmult._scripts.double_forest_polynomial import double_forest_polynomial, extract_double_forest_coefficients, extract_double_forest_coefficient
    # from schubmult.abc import x,y

    poly = double_forest_polynomial([2,0,1], lambda i: x[i], lambda j: y[j])
    # coeffs = extract_double_forest_coefficients(poly, x_genset=x, length=4)
    # print('coeffs:\n', coeffs)
    # print('main coeff:', coeffs.get((2,0,1,0)))
    # print('single helper:', extract_double_forest_coefficient(poly, [2,0,1], x_genset=x, length=4))
