from schubmult import abc
from schubmult.rings.abstract_schub_poly import GenericPrintingTerm
from schubmult.rings.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.symbolic import S, expand_seq


def monomial_quasisym(comp, length, genset):
    if any(c == 0 for c in comp):
        return S.Zero
    if len(comp) == length:
        return expand_seq(comp, genset)
    if len(comp) == 0:
        return S.One
    ret = S.Zero

    ret += monomial_quasisym(comp, length - 1, genset) + monomial_quasisym(comp[:-1], length - 1, genset) * genset[length] ** comp[-1]
    return ret


def stuffle(alpha, beta):
    """
    Computes the stuffle product of two compositions alpha and beta.
    Returns a dictionary where keys are resulting compositions (tuples)
    and values are their coefficients.
    """
    from collections import Counter

    # Base cases: if one is empty, the result is the other
    if not alpha:
        return {beta: 1}
    if not beta:
        return {alpha: 1}

    a, rest_a = alpha[0], alpha[1:]
    b, rest_b = beta[0], beta[1:]

    res = Counter()

    # 1. Take first from alpha: a + (rest_alpha * beta)
    for comp, coeff in stuffle(rest_a, beta).items():
        res[(a, *comp)] += coeff

    # 2. Take first from beta: b + (alpha * rest_beta)
    for comp, coeff in stuffle(alpha, rest_b).items():
        res[(b, *comp)] += coeff

    # 3. Stuff (sum) both: (a+b) + (rest_alpha * rest_beta)
    for comp, coeff in stuffle(rest_a, rest_b).items():
        res[(a + b, *comp)] += coeff

    return dict(res)


def quasi_schur_to_monomial(comp):
    """
    Computes the quasi-Schur function for composition comp in the monomial basis.
    Returns a dictionary where keys are compositions (tuples) and values are coefficients.

    Uses the standard composition tableau definition: sum over all descent compositions
    of standard composition tableaux of the given shape.
    """
    from itertools import permutations

    if not comp or all(c == 0 for c in comp):
        return {(): 1}

    # Total number of cells
    n = sum(comp)
    if n == 0:
        return {(): 1}

    # Build all standard composition tableaux of shape comp
    # A tableau is valid if entries increase along rows and down columns
    def descent_composition(tableau):
        """Extract the descent composition from a tableau (read row by row)."""
        flat = []
        for row in tableau:
            flat.extend(row)

        # Find descents in the reading word
        result = []
        current_run = 1
        for i in range(len(flat) - 1):
            if flat[i] > flat[i + 1]:
                result.append(current_run)
                current_run = 1
            else:
                current_run += 1
        result.append(current_run)
        return tuple(result)

    def is_valid_tableau(tableau):
        """Check if tableau has strictly increasing rows and weakly increasing columns."""
        # Check rows are strictly increasing
        for row in tableau:
            if any(row[i] >= row[i + 1] for i in range(len(row) - 1)):
                return False

        # Check columns are weakly increasing
        for col_idx in range(len(tableau[0]) if tableau else 0):
            for row_idx in range(len(tableau) - 1):
                if col_idx < len(tableau[row_idx]) and col_idx < len(tableau[row_idx + 1]):
                    if tableau[row_idx][col_idx] > tableau[row_idx + 1][col_idx]:
                        return False
        return True

    # Generate all tableaux by trying all ways to fill with 1..n
    from collections import Counter

    descent_comps = Counter()

    # Try all permutations and see which ones give valid tableaux
    for perm in permutations(range(1, n + 1)):
        # Build tableau from permutation
        tableau = []
        idx = 0
        for row_len in comp:
            tableau.append(list(perm[idx : idx + row_len]))
            idx += row_len

        if is_valid_tableau(tableau):
            desc_comp = descent_composition(tableau)
            descent_comps[desc_comp] += 1

    return dict(descent_comps)


class QSymElement(BaseSchubertElement):
    def expand(self, num_vars):
        """Expand the quasi-symmetric function in the given number of variables."""
        return sum([coeff * monomial_quasisym(comp, num_vars, self.ring.genset) for comp, coeff in self.items()])


class QSym(BaseSchubertRing):
    def __init__(self, genset=abc.x, domain=None):
        super().__init__(genset, None, domain)
        self.dtype = type("QSymElement", (QSymElement,), {"ring": self})

    def mul_pair(self, a, b):
        # a and b are compositions (tuples)
        return stuffle(a, b)

    def mul(self, a, b):
        # a and b are QSymElements
        result_monomials = self.zero
        for comp_a, coeff_a in a.items():
            for comp_b, coeff_b in b.items():
                for comp_c, coeff_c in self.mul_pair(comp_a, comp_b).items():
                    result_monomials += coeff_a * coeff_b * coeff_c * self(*comp_c)
        return result_monomials

    def expand(self, num_vars, *_, **__):
        return sum([coeff * monomial_quasisym(comp, num_vars, self.genset) for comp, coeff in self.items()])

    def rmul(self, a, b):
        return self.from_dict({comp: coeff * b for comp, coeff in a.items()})

    def printing_term(self, comp):
        return GenericPrintingTerm(tuple(comp), f"M{self.genset.label}")

    def new(self, *x):
        return self.from_dict({x: 1})

    def quasi_schur(self, *comp):
        """
        Returns the quasi-Schur function for the given composition
        expressed in the monomial basis.

        Args:
            *comp: A composition (tuple or sequence of positive integers)

        Returns:
            QSymElement representing the quasi-Schur function in monomial basis

        Example:
            >>> QS = QSym()
            >>> QS.quasi_schur(2, 1)  # quasi-Schur function for composition (2,1)
        """
        if len(comp) == 1 and isinstance(comp[0], (tuple, list)):
            comp = comp[0]

        monomial_expansion = quasi_schur_to_monomial(tuple(comp))
        return self.from_dict(monomial_expansion)

    # def from_expr(self, expr):
    #     dct = genset_dict_from_expr(expr, self.genset)
