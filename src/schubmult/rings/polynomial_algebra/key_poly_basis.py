from schubmult.combinatorics.permutation import Permutation
from schubmult.combinatorics.plactic import Plactic
from schubmult.symbolic import S, expand_seq
from schubmult.utils.perm_utils import add_perm_dict_with_coeff
from schubmult.utils.tuple_utils import pad_tuple

from ..printing import GenericPrintingTerm
from .base_polynomial_basis import PolynomialBasis

"""
Key polynomial (Demazure character) basis for Schubert calculus.

This module implements the key polynomial basis, computed via Demazure
operators on plactic (Yamanouchi) tableaux.
"""

def traverse_demaz(pl, w):
    """Traverse the Demazure graph starting from plactic element *pl* along the code word of *w*.

    Yields all distinct elements reachable by successive lowering operators.
    """
    stack = [(pl, 0)]
    word = list(reversed(w.code_word))
    while stack:
        current, index = stack.pop()
        yield current
        for index2 in range(index, len(word)):
            desc2 = word[index2]
            next_elem = current.lowering_operator(desc2)
            if next_elem is not None:
                stack.append((next_elem, index2))




def _key_polynomial(comp, genset):
    """Compute the key polynomial (Demazure character) for a weak composition."""
    w = Permutation.sorting_perm(comp, reverse=True)
    hw = sorted(comp, reverse=True)
    top_pl = Plactic.yamanouchi(hw)
    visited = set()
    ret = S.Zero
    for pl in traverse_demaz(top_pl, w):
        if pl in visited:
            continue
        visited.add(pl)
        ret += expand_seq(pl.crystal_weight, genset)
    return ret

class KeyPolyBasis(PolynomialBasis):
    """Key polynomial (Demazure character) basis.

    Keys are weak compositions. Key polynomials are characters of Demazure
    modules, computed by applying Demazure operators to highest-weight
    plactic tableaux.
    """
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"KeyPoly{k}", "")

    @classmethod
    def dual_basis(cls):
        """Return the dual free algebra basis class (:class:`KeyBasis`)."""
        from ..free_algebra.key_basis import KeyBasis
        return KeyBasis

    def __init__(self, genset):
        from .monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        """Expand a key polynomial key into a dict of monomial exponent tuples."""
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_key_polynomial(key, self.genset), self.genset).items()}
        return dct

    def expand(self, dct):
        """Expand a key basis dict into a symbolic polynomial expression."""
        return sum([v * _key_polynomial(key, self.genset) for key, v in dct.items()])

    def transition_monomial(self, dct):
        """Transition from key basis to monomial basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition(self, other_basis):
        """Return a transition function from key basis to *other_basis*."""
        from .monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    @property
    def zero_monom(self):
        return self.as_key([])
