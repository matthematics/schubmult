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

    ret += monomial_quasisym(comp, length - 1, genset) + monomial_quasisym(comp[:-1], length - 1, genset) * genset[length]**comp[-1]
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


class QSymElement(BaseSchubertElement):
    pass


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

    # def from_expr(self, expr):
    #     dct = genset_dict_from_expr(expr, self.genset)
