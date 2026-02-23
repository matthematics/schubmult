from schubmult.combinatorics.indexed_forests import weak_composition_to_indfor
from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict_with_coeff


def _fundamental_slide_polynomial(comp, genset):
    compat_seq = []
    for i, c in enumerate(comp):
        compat_seq.extend([i + 1] * c)
    compat_seq.reverse()
    return _compat_seq_poly(tuple(compat_seq), genset)


def decreasing_labelings(root, max_val, used_vals=None):
    def _subtree_size(node):
        if node is None:
            return 0
        return 1 + _subtree_size(node.left) + _subtree_size(node.right)

    ret = set()
    used = set() if used_vals is None else set(used_vals)
    required_size = _subtree_size(root)
    available_labels = [val for val in range(1, max_val + 1) if val not in used]
    if len(available_labels) < required_size:
        return ret

    if root.left is None and root.right is None:
        return {(val,) for val in available_labels}

    left_size = _subtree_size(root.left)
    right_size = _subtree_size(root.right)

    for val in available_labels:
        remaining_smaller = sum(1 for a in available_labels if a < val)
        if remaining_smaller < left_size + right_size:
            continue

        vals_used = used | {val}
        if root.left is None:
            right_labelings = decreasing_labelings(root.right, val - 1, used_vals=vals_used)
            for right in right_labelings:
                ret.add((val, *right))
            continue

        if root.right is None:
            left_labelings = decreasing_labelings(root.left, val - 1, used_vals=vals_used)
            for left in left_labelings:
                ret.add((*left, val))
            continue

        right_labelings = decreasing_labelings(root.right, val - 1, used_vals=vals_used)
        for right in right_labelings:
            new_vals_used = vals_used | set(right)
            left_labelings = decreasing_labelings(root.left, val - 1, used_vals=new_vals_used)
            for left in left_labelings:
                ret.add((*left, val, *right))
    return ret
    # right_labelings = set()
    # left_labelings = set()

    # if root.left is not None:
    #     left_labelings.update(decreasing_labelings(root.left, val - 1))
    # if root.right is None:
    #     ret.update(tuple(list(left) + [val]) for left in left_labelings)
    #     return ret
    # if root.left is None:
    #     ret.update(tuple([val] + list(right)) for right in right_labelings)
    #     return ret
    # ret.update(tuple(list(left) + [val] + list(right)) for left in left_labelings for right in right_labelings)
    # return ret


def word_from_labeling(root, labeling):
    from schubmult import Permutation

    trav = list(root.inorder_traversal)
    assert len(trav) == len(labeling)
    perm = Permutation(labeling)
    # print(f"Labeling: {labeling}, Perm: {perm}, Rhos: {[node.rho for node in trav]}")
    return tuple(trav[(~perm)[i] - 1].rho for i in range(len(trav)))


def word_to_canonical_composition_accurate(word):
    """
    Converts a word to its canonical weak composition by first
    reducing it to its unique non-increasing 'slide' equivalent.
    """
    if not word:
        return []

    # Slide logic: We can decrement an entry i_j to i_j - 1
    # if it doesn't create a new descent or violate a_j <= i_j.
    # For 423 -> 422, the '3' can slide down to '2' because
    # it is preceded by a '2' and 2 <= 3.
    canonical_word = list(word)
    done_any = True
    while done_any:
        done_any = False
        for i in range(len(canonical_word) - 1, 0, -1):
            # If a number is greater than the one to its left,
            # it can often be 'slid' down.
            if canonical_word[i - 1] < canonical_word[i]:
                canonical_word[i] = canonical_word[i - 1]
                done_any = True
    # check same descents as word
    for i in range(len(word) - 1):
        if (word[i] > word[i + 1]) != (canonical_word[i] > canonical_word[i + 1]):
            return None
    # Sort non-increasingly to get the standard representative W_c
    # canonical_word.sort(reverse=True)
    # Now convert the canonical word to a composition
    max_val = max(canonical_word)
    composition = [0] * max_val
    for val in canonical_word:
        composition[val - 1] += 1
    return tuple(composition)


def _forest_polynomial_from_indfor(indfor, genset):
    if len(indfor) == 0:
        return S.One

    ret = S.One

    for root in indfor:
        ret *= _forest_polynomial_from_root(root, genset)
    return ret


def _forest_polynomial_from_root(root, genset, min_value=1):
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

def _flatten_seq_to_comp(seq, length):
    if len(seq) == 0:
        return ()
    comp_seq = [*seq[:1]]
    for i in range(1, len(seq)):
        if seq[i] < seq[i - 1]:
            comp_seq.append(min(seq[i], comp_seq[i - 1] - 1))
        #if seq[i] > seq[i - 1]:
        else:
            comp_seq.append(comp_seq[i - 1])
    if not all(a > 0 for a in comp_seq):
        return None
    composition = [0] * length
    for val in comp_seq:
        composition[val - 1] += 1
    return tuple(composition)
class ForestPolyBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list)

    def as_key(self, x):
        return tuple(x)

    def printing_term(self, k):
        return GenericPrintingTerm(f"Forest{k}", "")

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
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        self._genset = genset
        self._monomial_basis = MonomialBasis(genset=self.genset)

    def to_monoms(self, key):
        from schubmult.symbolic.poly.variables import genset_dict_from_expr

        def pad_tuple(tup, length):
            return (*tup, *(0,) * (length - len(tup)))

        dct = {pad_tuple(k, len(key)): v for k, v in genset_dict_from_expr(_forest_polynomial_from_indfor(weak_composition_to_indfor(key), self.genset), self.genset).items()}
        return dct

    def expand(self, dct):
        return sum([v * _forest_polynomial_from_indfor(weak_composition_to_indfor(key), self.genset) for key, v in dct.items()])

    def transition_monomial(self, dct):
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_monoms(k), coeff=v)
        return res

    def transition_fundamental_slide(self, dct):
        res = {}
        for k, v in dct.items():
            res = add_perm_dict_with_coeff(res, self.to_fundamental_slide(k), coeff=v)
        return res

    def to_fundamental_slide(self, key):
        indfor = weak_composition_to_indfor(key)
        if len(indfor) > 1:
            raise NotImplementedError("Transition to fundamental slide for products not implemented yet")
        dct = {}
        for root in indfor:
            size = len(list(root.inorder_traversal))
            for labeling in decreasing_labelings(root, size):
                word = word_from_labeling(root, labeling)
                comp = _flatten_seq_to_comp(word, len(key))
                if comp is None:
                    continue
                # print(f"Word: {word}, Comp: {comp}, Forest Poly: {_forest_polynomial_from_root(root, self.genset)}, FSlide: {_fundamental_slide_polynomial(comp, self.genset)}, Compat Seq Poly: {_compat_seq_poly(word, self.genset)}")
                assert (_fundamental_slide_polynomial(comp, self.genset) - _compat_seq_poly(word, self.genset)).expand() == 0, f"Word: {word}, Comp: {comp}, Forest Poly: {_forest_polynomial_from_root(root, self.genset)}, FSlide: {_fundamental_slide_polynomial(comp, self.genset)}, Compat Seq Poly: {_compat_seq_poly(word, self.genset)}"
                dct[comp] = dct.get(comp, S.Zero) + S.One
        return dct

    def transition(self, other_basis):
        from schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis import FundamentalSlidePolyBasis
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial
        if isinstance(other_basis, FundamentalSlidePolyBasis):
            return self.transition_fundamental_slide

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def from_expr(self, expr):
        dct = self.monomial_basis.from_expr(expr)
        return self.monomial_basis.transition_slide(self, dct)

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
