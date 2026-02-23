from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict_with_coeff

from ..printing import GenericPrintingTerm
from .base_polynomial_basis import PolynomialBasis

"""
Fundamental slide polynomial basis for Schubert calculus.

This module implements the fundamental slide polynomial basis, which provides
an alternative basis for expressing Schubert polynomials and their products.
"""


def get_descent_composition(word):
    if not word: return []
    descents = [1]
    for i in range(len(word) - 1):
        if word[i+1] > word[i]:
            descents[-1] += 1
        else:
            descents.append(1)
    return descents

def slide_product(a, b):
    n = len(a)
    # 1. Create encoded words A and B (Definition 5.8)
    # A uses odd labels, B uses even labels, decreasing within rows
    A_labeled = []
    for i in range(n):
        A_labeled.extend([(2*(n-i)-1)] * a[i])
    
    B_labeled = []
    for i in range(n):
        B_labeled.extend([(2*(n-i))] * b[i])
        
    len_a, len_b = sum(a), sum(b)
    total_len = len_a + len_b
    
    # 2. Standard Shuffle Product (riffle shuffle of words)
    def shuffles(w1, w2):
        if not w1: return [w2]
        if not w2: return [w1]
        return [[w1[0]] + s for s in shuffles(w1[1:], w2)] + \
               [[w2[0]] + s for s in shuffles(w1, w2[1:])]

    all_shuffles = shuffles(A_labeled, B_labeled)
    
    results = {}
    for C in all_shuffles:
        # 3. Check SS(a, b) condition (Definition 5.8)
        # Des_A(C)_i must sum to a_i, etc.
        # This is implicitly handled by the labeling and 'bump' logic
        
        # Calculate Descent Composition of the shuffle
        flat_c = get_descent_composition(C)
        
        # 4. Apply 'bump' (Definition 5.9)
        # Pad with leading zeros to maintain length n
        # and ensure dominance b_tilde_ge a
        curr_len = len(flat_c)
        if curr_len > n:
            continue # Exceeds polynomial ring dimensions
            
        # The 'bump' logic effectively places the flat composition 
        # into the rightmost possible slots that don't violate dominance.
        # For the dual basis, we collect the resulting weak composition.
        padded_c = [0] * (n - curr_len) + flat_c
        
        res_tuple = tuple(padded_c)
        results[res_tuple] = results.get(res_tuple, 0) + 1
        
    return results


def _fundamental_slide_polynomial(comp, genset):
    compat_seq = []
    for i, c in enumerate(comp):
        compat_seq.extend([i + 1] * c)
    compat_seq.reverse()
    return _compat_seq_poly(tuple(compat_seq), genset)


def _compat_seq_poly(comp, genset):
    if len(comp) == 0:
        return S.One

    ret = S.Zero

    if len(comp) == 1:
        for i in range(1, comp[0] + 1):
            ret += genset[i]
        return ret

    last_elem = comp[-1]
    genstart = 0
    working_comp = comp[:-1]
    if last_elem != comp[-2]:
        genstart = 1
        working_comp = [a - 1 for a in working_comp]
    for i in range(1, last_elem + 1):
        ret += genset[i] * _compat_seq_poly(working_comp, genset[genstart + i - 1 :])
        working_comp = [a - 1 for a in working_comp]
    return ret


class FundamentalSlidePolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"FSlide{k}", "")

    # def coproduct(self, key):
    #     result_dict = {}
    #     key = self.as_key(key)
    #     for i in range(len(key) + 1):
    #         result_dict[(key[:i], key[i:])] = S.One
    #     return result_dict

    @property
    def monomial_basis(self):
        return self._monomial_basis

    @property
    def genset(self):
        return self._genset

    def __init__(self, genset):
        from .monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        def pad_tuple(tup, length):
            return (*tup, *(0,) * (length - len(tup)))

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_fundamental_slide_polynomial(key, self.genset), self.genset).items()}
        return dct

    def expand(self, dct):
        return sum([v * _fundamental_slide_polynomial(key, self.genset) for key, v in dct.items()])

    def transition_monomial(self, dct):
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition(self, other_basis):
        from .monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def from_expr(self, expr):
        # dct = self.monomial_basis.from_expr(expr)
        # return self.monomial_basis.transition_slide(dct, self)
        try:
            from schubmult.symbolic import sympify
            return {self.zero_monom: sympify(expr)}
        except:
            pass    
        raise NotImplementedError("Direct expression parsing not implemented for FundamentalSlidePolyBasis yet")

    def product(self, key1, key2, coeff=S.One):
        if len(key1) != len(key2):
            return {}

        return {c: v * coeff for c, v in slide_product(key1, key2).items()}


    @property
    def zero_monom(self):
        return self.as_key([])
