from schubmult.rings.polynomial_algebra.base_polynomial_basis import PolynomialBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict_with_coeff

"""
Fundamental slide polynomial basis for Schubert calculus.

This module implements the fundamental slide polynomial basis, which provides
an alternative basis for expressing Schubert polynomials and their products.
"""


def get_descent_composition(word):
    if not word:
        return []
    descents = [1]
    for i in range(len(word) - 1):
        if word[i+1] > word[i]:
            descents[-1] += 1
        else:
            descents.append(1)
    return descents

def slide_product(a, b):
    n = max(len(a), len(b))
    if len(a) < n:
        a = (*a, *(0 for _ in range(n - len(a))))
    if len(b) < n:
        b = (*b, *(0 for _ in range(n - len(b))))

    def _prefix_dominates(x, y):
        if sum(x) != sum(y):
            return False
        sx = 0
        sy = 0
        for i in range(max(len(x), len(y))):
            sx += x[i] if i < len(x) else 0
            sy += y[i] if i < len(y) else 0
            if sx < sy:
                return False
        return True

    def _run_start_indices(word_vals):
        if not word_vals:
            return []
        starts = [0]
        for i in range(1, len(word_vals)):
            if word_vals[i - 1] > word_vals[i]:
                starts.append(i)
        return starts

    def _run_stats(tagged_word):
        vals = [v for v, _ in tagged_word]
        starts = _run_start_indices(vals)
        if not starts:
            return (), (), ()
        starts.append(len(tagged_word))
        des_total = []
        des_a = []
        des_b = []
        for i in range(len(starts) - 1):
            lo = starts[i]
            hi = starts[i + 1]
            run = tagged_word[lo:hi]
            des_total.append(len(run))
            des_a.append(sum(1 for _, tag in run if tag == "A"))
            des_b.append(sum(1 for _, tag in run if tag == "B"))
        return tuple(des_total), tuple(des_a), tuple(des_b)

    def _insert_zeros_into_comp(comp, zero_count, length):
        if len(comp) + zero_count != length:
            return
        if not comp:
            yield tuple(0 for _ in range(length)), ()
            return
        from itertools import combinations

        for nonzero_positions in combinations(range(length), len(comp)):
            out = [0] * length
            for i, pos in enumerate(nonzero_positions):
                out[pos] = comp[i]
            yield tuple(out), tuple(nonzero_positions)

    def _place_by_positions(comp, positions, length):
        out = [0] * length
        for i, pos in enumerate(positions):
            out[pos] = comp[i]
        return tuple(out)

    def _shuffle_tagged_words(w1, w2):
        if not w1:
            yield tuple(w2)
            return
        if not w2:
            yield tuple(w1)
            return
        for rest in _shuffle_tagged_words(w1[1:], w2):
            yield (w1[0], *rest)
        for rest in _shuffle_tagged_words(w1, w2[1:]):
            yield (w2[0], *rest)

    def _strictly_dominates(x, y):
        return _prefix_dominates(x, y) and not _prefix_dominates(y, x)

    def _bump(des_total, des_a, des_b):
        zero_count = n - len(des_total)
        if zero_count < 0:
            return None

        best = None
        for candidate, positions in _insert_zeros_into_comp(des_total, zero_count, n):
            candidate_a = _place_by_positions(des_a, positions, n)
            candidate_b = _place_by_positions(des_b, positions, n)
            if not (_prefix_dominates(candidate_a, a) and _prefix_dominates(candidate_b, b)):
                continue
            if best is None:
                best = candidate
                continue
            if _strictly_dominates(best, candidate):
                best = candidate
            elif (not _prefix_dominates(candidate, best)) and (not _prefix_dominates(best, candidate)):
                if candidate < best:
                    best = candidate
        return best

    a_word = []
    for i in range(n):
        a_word.extend([(2 * (n - i) - 1, "A")] * a[i])

    b_word = []
    for i in range(n):
        b_word.extend([(2 * (n - i), "B")] * b[i])

    results = {}
    for c_word in _shuffle_tagged_words(tuple(a_word), tuple(b_word)):
        des_total, des_a, des_b = _run_stats(c_word)
        if not (_prefix_dominates(des_a, a) and _prefix_dominates(des_b, b)):
            continue
        bumped = _bump(des_total, des_a, des_b)
        if bumped is None:
            continue
        results[bumped] = results.get(bumped, 0) + 1

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
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

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
        from schubmult.rings.polynomial_algebra.monomial_basis import MonomialBasis

        if isinstance(other_basis, MonomialBasis):
            return self.transition_monomial

        return lambda x: PolynomialBasis.compose_transition(self.monomial_basis.transition(other_basis), self.transition_monomial(x))

    def from_expr(self, expr):
        # dct = self.monomial_basis.from_expr(expr)
        # return self.monomial_basis.transition_slide(dct, self)
        try:
            from schubmult.symbolic import sympify
            return {self.zero_monom: sympify(expr)}
        except Exception:
            pass
        raise NotImplementedError("Direct expression parsing not implemented for FundamentalSlidePolyBasis yet")

    def product(self, key1, key2, coeff=S.One):
        return {c: v * coeff for c, v in slide_product(key1, key2).items()}


    @property
    def zero_monom(self):
        return self.as_key([])


if __name__ == "__main__":
    from schubmult import PolynomialAlgebra, Sx

    FSlide = PolynomialAlgebra(FundamentalSlidePolyBasis(Sx.genset))
    print(FSlide(0, 2, 0, 3) * FSlide(1, 0, 0, 1))
