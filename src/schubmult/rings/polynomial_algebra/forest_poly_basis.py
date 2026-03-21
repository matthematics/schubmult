from schubmult.combinatorics.indexed_forests import decreasing_labelings, weak_composition_to_indfor, word_from_labeling
from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict_with_coeff
from schubmult.utils.tuple_utils import pad_tuple


def _forest_polynomial_from_indfor(indfor, genset):
    """Compute the forest polynomial for an indexed forest by multiplying root contributions."""
    if len(indfor) == 0:
        return S.One

    ret = S.One

    for root in indfor:
        ret *= _forest_polynomial_from_root(root, genset)
    return ret


def _forest_polynomial_from_root(root, genset, min_value=1):
    """Recursively compute the polynomial contribution from a single indexed forest root."""
    if root.rho == 0:
        return S.One
    if root.left is None and root.right is None:
        return sum([genset[i] for i in range(min_value, root.rho + 1)])
    ret = S.Zero
    for val in range(min_value, root.rho + 1):
        term = genset[val]
        if root.left is not None:
            term *= _forest_polynomial_from_root(root.left, genset, val)
        if root.right is not None:
            term *= _forest_polynomial_from_root(root.right, genset, val + 1)
        ret += term
    return ret


def _flatten_seq_to_comp(seq, length):
    """Convert a value sequence to a weak composition by run-length encoding with dominance."""
    if len(seq) == 0:
        return ()
    comp_seq = [*seq[:1]]
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            comp_seq.append(min(seq[i], comp_seq[i - 1] - 1))
        # if seq[i] > seq[i - 1]:
        else:
            comp_seq.append(comp_seq[i - 1])
    if not all(a > 0 for a in comp_seq):
        return None
    composition = [0] * length
    for val in comp_seq:
        composition[val - 1] += 1
    return tuple(composition)


class ForestPolyBasis(PolynomialBasis):
    """Forest polynomial basis.

    Keys are weak compositions encoding indexed forests. Forest polynomials
    are computed by summing over decreasing labelings of the corresponding
    forest structure.
    """
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"Forest{k}", "")

    def __init__(self, genset):
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis
        super().__init__(genset=genset)
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        """Expand a forest key into a dict of monomial exponent tuples."""
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_forest_polynomial_from_indfor(weak_composition_to_indfor(key), self.genset), self.genset).items()}
        return dct

    def expand(self, dct):
        """Expand a forest basis dict into a symbolic polynomial expression."""
        return sum([v * _forest_polynomial_from_indfor(weak_composition_to_indfor(key), self.genset) for key, v in dct.items()])

    def transition_monomial(self, dct):
        """Transition from forest basis to monomial basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition_fundamental_slide(self, dct):
        """Transition from forest basis to fundamental slide basis."""
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_fundamental_slide(k), coeff=v)
        return res

    def to_fundamental_slide(self, key):
        """Express a single forest key in the fundamental slide basis."""
        from schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis import slide_product
        indfor = weak_composition_to_indfor(key)

        ret = None
        if len(indfor) == 0:
            return {key: S.One}
        for root in indfor:
            dct = {}
            size = len(list(root.inorder_traversal))
            for labeling in decreasing_labelings(root, size):
                word = word_from_labeling(root, labeling)
                comp = _flatten_seq_to_comp(word, len(key))
                if comp is None:
                    continue
                dct[comp] = dct.get(comp, S.Zero) + S.One
            if ret is None:
                ret = dct
            else:
                new_ret = {}
                for k1, v1 in ret.items():
                    for k2, v2 in dct.items():
                        s_prod = slide_product(k1, k2)
                        for comp, coeff in s_prod.items():
                            new_ret[comp] = new_ret.get(comp, S.Zero) + v1 * v2 * coeff
                ret = new_ret
        return ret

    @classmethod
    def dual_basis(cls):
        """Return the dual free algebra basis class (:class:`ForestBasis`)."""
        from ..free_algebra.forest_basis import ForestBasis
        return ForestBasis

    def transition(self, other_basis):
        """Return a transition function from forest basis to *other_basis*."""
        from schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis import FundamentalSlidePolyBasis
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial
        if isinstance(other_basis, FundamentalSlidePolyBasis):
            return self.transition_fundamental_slide

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def product(self, key1, key2, coeff=S.One):
        """Multiply two forest keys by transitioning through the Schubert basis."""
        from schubmult.utils.perm_utils import add_perm_dict

        from ..schubert.schubert_ring import Sx
        from .schubert_poly_basis import SchubertPolyBasis

        schub = SchubertPolyBasis(Sx)
        left = self.transition(schub)({key1: coeff})
        right = self.transition(schub)({key2: S.One})

        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, schub.transition(self)(schub.product(key_schub_left, key_schub_right, v * v2)))
        return ret

    @property
    def zero_monom(self):
        return self.as_key([])


if __name__ == "__main__":
    from symengine import expand

    from schubmult import PolynomialAlgebra, Sx
    from schubmult.rings.polynomial_algebra import FundamentalSlidePolyBasis

    For = PolynomialAlgebra(ForestPolyBasis(Sx.genset))
    print("Fundy ", For([0, 2, 0, 1]).change_basis(FundamentalSlidePolyBasis(Sx.genset)).expand().expand())
    print("Monny, ", For([0, 2, 0, 1]).expand().expand())
    print(expand(For([0, 2, 0, 1]).change_basis(FundamentalSlidePolyBasis(Sx.genset)).expand() - For([0, 2, 0, 1]).expand()))
