from functools import cache

import schubmult.rings.free_algebra as fa
from schubmult.combinatorics.permutation import uncode
from schubmult.symbolic import S, Symbol
from schubmult.utils.perm_utils import add_perm_dict
from schubmult.utils.tuple_utils import pad_tuple

from ..schubert.schubert_ring import DSx, Sx
from ..schubert.separated_descents import SeparatedDescentsRing
from .free_algebra_basis import FreeAlgebraBasis

splugSx = SeparatedDescentsRing(Sx([]).ring)
ADSx = SeparatedDescentsRing(DSx([]).ring)


class WordBasis(FreeAlgebraBasis):
    """Word basis of the free algebra.

    Keys are tuples of nonnegative integers representing words. This is the
    fundamental basis through which all other bases perform their operations
    via basis transitions.
    """

    @staticmethod
    def _weak_compositions(length, total):
        """Yield all weak compositions of *total* into *length* nonneg parts."""
        if length == 0:
            if total == 0:
                yield ()
            return
        if length == 1:
            yield (total,)
            return
        for i in range(total + 1):
            for tail in WordBasis._weak_compositions(length - 1, total - i):
                yield (i, *tail)

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a tuple or list."""
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* to a tuple key."""
        return tuple(x)

    @classmethod
    def from_rc_graph(cls, rc_graph):
        """Return the length vector of the RC graph as a word key."""
        return {rc_graph.length_vector(): 1}

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        """Concatenate two words."""
        return {(*key1, *key2): coeff}

    @classmethod
    def inject(cls, key1, i, key2, coeff=S.One):
        """Insert *key2* into *key1* at position *i*."""
        if i > len(key1):
            raise IndexError(f"Insertion index {i} out of range for word of length {len(key1)}")
        return {(*key1[:i], *key2, *key1[i:]): coeff}

    @classmethod
    def prefix(cls, key, length, coeff=S.One):
        """Return the first *length* letters of *key*."""
        return {key[:length]: coeff}

    @classmethod
    def suffix(cls, key, length, coeff=S.One):
        """Return the last *length* letters of *key*."""
        return {key[len(key) - length:]: coeff}

    @classmethod
    def interval(cls, key, start, stop, coeff=S.One):
        """Return the subword ``key[start:stop]``."""
        return {key[start:stop]: coeff}

    zero_monom = ()

    @classmethod
    @cache
    def coproduct(cls, key, coeff=S.One):
        """Compute the additive coproduct of a word.

        Decomposes each letter into all (i, key-i) splittings and combines
        via a divide-and-conquer tensor product.
        """
        if len(key) == 0:
            return {((), ()): coeff}
        if len(key) == 1:
            key = key[0]
            dct = {}
            for i in range(key + 1):
                dct[((i,), (key - i,))] = coeff
            return dct
        mid = len(key) // 2
        cp1 = cls.coproduct(key[:mid], coeff)
        cp2 = cls.coproduct(key[mid:])
        ret = {}
        for k0, v0 in cp1.items():
            for k1, v1 in cp2.items():
                ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
        return ret

    @classmethod
    @cache
    def bcoproduct(cls, key, coeff=S.One):
        """Compute the bar-coproduct of a word.

        Like :meth:`coproduct` but drops empty factors (zeros map to
        the empty tuple).
        """
        if len(key) == 0:
            return {((), ()): coeff}
        if len(key) == 1:
            key = key[0]
            dct = {}
            for i in range(key + 1):
                dct[((i,) if i != 0 else (), (key - i,) if key - i != 0 else ())] = coeff
            return dct
        mid = len(key) // 2
        cp1 = cls.bcoproduct(key[:mid], coeff)
        cp2 = cls.bcoproduct(key[mid:])
        ret = {}
        for k0, v0 in cp1.items():
            for k1, v1 in cp2.items():
                ret = add_perm_dict(ret, {((*k0[0], *k1[0]), (*k0[1], *k1[1])): v0 * v1})
        return ret

    @classmethod
    def try_internal_product(cls, key1, key2, coeff=S.One):
        """Compute the internal product via integer matrices (requires SageMath).

        Uses shifted keys (incremented by 1) with ``IntegerMatrices``.
        """
        from sage.combinat.integer_matrices import IntegerMatrices

        bkey1 = [a + 1 for a in key1]
        bkey2 = [a + 1 for a in key2]
        mats = IntegerMatrices(bkey1, bkey2).list()

        def mat_to_key(mat):
            res = []
            for row in mat:
                for val in row:
                    if val != 0:
                        res.append(val - 1)
            return tuple(res)

        ret_dict = {}
        for mat in mats:
            key = mat_to_key(mat)
            if key not in ret_dict:
                ret_dict[key] = S.Zero
            ret_dict[key] += coeff
        return ret_dict

    @classmethod
    def internal_product(cls, key1, key2, coeff=S.One):
        """Compute the internal product of two words via integer matrices (requires SageMath).

        Words must not contain zeros. Returns the dict of result
        words weighted by *coeff*.
        """
        if 0 in key1 or 0 in key2:
            return {}
        from sage.combinat.integer_matrices import IntegerMatrices

        mats = list(IntegerMatrices(key1, key2))

        def mat_to_key(mat):
            res = []
            for row in mat:
                for val in row:
                    if val != 0:
                        res.append(int(val))
            return tuple(res)

        ret_dict = {}
        for mat in mats:
            key = mat_to_key(mat)
            if key not in ret_dict:
                ret_dict[key] = S.Zero
            ret_dict[key] += coeff
        return ret_dict

    @classmethod
    def printing_term(cls, k):
        """Return a bracket-notation symbol like ``[210]`` for the word *k*."""
        if all(a < 10 for a in k):
            return Symbol("[" + "".join([str(a) for a in k]) + "]")
        return Symbol("[" + " ".join([str(a) for a in k]) + "]")

    @staticmethod
    @cache
    def tup_expand(tup):
        """Expand a word tuple into the single Schubert basis via divide-and-conquer."""
        res = splugSx([])
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            return splugSx(uncode(tup), 1)
        mid = len(tup) // 2
        return WordBasis.tup_expand(tup[:mid]) * WordBasis.tup_expand(tup[mid:])

    @staticmethod
    @cache
    def jbasis_tup_expand(tup):
        """Expand a word tuple into the Z basis."""
        from .z_basis import ZBasis

        JB = fa.FreeAlgebra(basis=ZBasis)
        res = JB()
        if len(tup) == 0:
            return res
        if len(tup) == 1:
            if tup[0] == 0:
                return res
            return JB(*tup)
        mid = len(tup) // 2
        return WordBasis.jbasis_tup_expand(tup[:mid]) * WordBasis.jbasis_tup_expand(tup[mid:])

    @classmethod
    def transition_schubert(cls, key):
        """Transition a word key to the Schubert basis."""
        return dict(WordBasis.tup_expand(key))

    @classmethod
    def transition_jbasis(cls, key):
        """Transition a word key to the J basis via Pieri products."""
        from .j_basis import JBasis
        from .schubert_basis import SchubertBasis

        t = S.One
        JB = fa.FreeAlgebra(basis=JBasis)
        ASx = fa.FreeAlgebra(basis=SchubertBasis)
        assert all(a != 0 for a in key), "Transition from WordBasis to JBasis is only implemented for keys with no zeros."
        pangea = JB()

        for a in reversed(key):
            new_pangea = JB.from_dict({})
            for bucket, coeff in pangea.items():
                the_pieri = coeff * ASx(uncode([a])) * ASx(uncode(bucket))
                for (perm, n), v in the_pieri.items():
                    if 0 in perm.trimcode:
                        new_pangea += v * t * JB(*[a for a in perm.trimcode if a != 0])
                    else:
                        new_pangea += v * JB(*perm.trimcode)
            pangea = new_pangea

        return pangea

    @classmethod
    def transition_jtbasis(cls, key):
        """Transition a word key to the JT basis via normalization."""
        from .jt_basis import JTBasis

        FA = fa.FreeAlgebra(WordBasis)
        JB = fa.FreeAlgebra(basis=JTBasis)
        res = FA.from_dict(JTBasis.normalize_dct(FA(*key)))

        ret = JB.from_dict({})

        while res != FA.zero:
            res = FA.from_dict({k: v for k, v in res.items() if v != S.Zero})
            tup = next(iter(sorted(res.keys())))
            c = res[tup]
            new_tup = tuple([int(a) for a in tup if a != 0])
            pw = len(tup) - len(new_tup)

            val = c * JB(new_tup, pw).change_basis(WordBasis)
            res -= val
            ret += c * JB(new_tup, pw)
        return dict(ret)

    @classmethod
    def transition_forest(cls, key):
        """Transition a word key to the forest basis via RC graph enumeration."""
        from schubmult.rings.combinatorial import RCGraphRing
        r = RCGraphRing()
        dct = {}
        all_rcs = r.monomial(*key)
        seen = {}
        for rc in all_rcs:
            indfor = rc.forest_invariant
            code_key = pad_tuple(indfor.forest.code, len(key))
            if code_key not in seen:
                seen[code_key] = indfor
            elif seen[code_key] != indfor:
                continue
            dct[code_key] = dct.get(code_key, S.Zero) + S.One
        return dct

    @classmethod
    def dual_basis(cls):
        """Return the MonomialBasis as the dual of WordBasis."""
        from ..polynomial_algebra.monomial_basis import MonomialBasis
        return MonomialBasis
    # @classmethod
    # def transition_forest(cls, key):
    #     from schubmult.combinatorics.indexed_forests import word_to_indexed_forest
    #     from schubmult.rings.combinatorial import RCGraphRing
        # r = RCGraphRing()
        # dct = {}
        # all_rcs = r.monomial(*key)

        # perms_used = set()

        # for rc in all_rcs:
        #     word = list(reversed(rc.perm_word))
        #     indfor = word_to_indexed_forest(word)
        #     assert len(indfor.code) == len(key), f"Length of indexed forest code {indfor.code} does not match length of key {key}."
        #     if rc.perm not in perms_used:
        #         dct[indfor.code] = dct.get(indfor.code, S.Zero) + S.One
        #         perms_used.add(rc.perm)
        # return dct

        # Pieri

    # @classmethod
    # def transition_fundamental_slide(cls, key):
    #     raise NotImplementedError("Transition from WordBasis to FundamentalSlideBasis is not implemented yet.")

    @classmethod
    def transition_monomial_slide(cls, key):
        """Transition a word key to the monomial slide basis."""
        from schubmult.abc import x
        from schubmult.rings.polynomial_algebra.monomial_slide_poly_basis import MonomialSlidePolyBasis

        key = tuple(key)
        length = len(key)
        total = sum(key)
        poly_basis = MonomialSlidePolyBasis(x)

        ret = {}
        for candidate in cls._weak_compositions(length, total):
            coeff = poly_basis.to_monoms(candidate).get(key, S.Zero)
            if coeff != S.Zero:
                ret[candidate] = coeff
        return ret

        # key = tuple(key)
        # n = len(key)
        # flat_key = tuple(a for a in key if a > 0)

        # if len(flat_key) == 0:
        #     return {tuple([0] * n): S.One}

        # prefix_key = []
        # running = 0
        # for a in key:
        #     running += a
        #     prefix_key.append(running)

        # ret = {}
        # m = len(flat_key)
        # for positions in combinations(range(n), m):
        #     candidate = [0] * n
        #     for i, pos in enumerate(positions):
        #         candidate[pos] = flat_key[i]

        #     running = 0
        #     dominates = True
        #     for i, a in enumerate(candidate):
        #         running += a
        #         if running < prefix_key[i]:
        #             dominates = False
        #             break

        #     if dominates:
        #         key_candidate = tuple(candidate)
        #         ret[key_candidate] = ret.get(key_candidate, S.Zero) + S.One

        # return ret

    @classmethod
    def transition_zbasis(cls, key):
        from .z_basis import ZBasis

        FA = fa.FreeAlgebra(WordBasis)
        JB = fa.FreeAlgebra(basis=ZBasis)
        res = FA(*key)

        ret = JB.from_dict({})
        while res != FA.zero:
            tup = next(iter(sorted(res.keys())))
            c = res[tup]
            val = c * JB(*tup).change_basis(WordBasis)
            res -= val
            ret += c * JB(*tup)
        return dict(ret)

    @classmethod
    def transition_nelementary(cls, tup):
        from sage.all import Composition

        pain = Composition(tup)
        piss = list(pain.finer())
        return {tuple([int(p) for p in pi]): (S.NegativeOne ** (sum(tup) - len(pi))) for pi in piss}

    @classmethod
    def transition_key(cls, key):
        from schubmult.rings.combinatorial import RCGraphRing
        r = RCGraphRing()
        dct = {}
        all_rcs = r.monomial(*key)
        highest_weights = set({rc.to_highest_weight()[0] for rc in all_rcs})
        seen = {}
        for hw in highest_weights:
            code_key = hw.extremal_weight
            if code_key not in seen:
                seen[code_key] = hw
            elif seen[code_key] != hw:
                continue
            dct[code_key] = dct.get(code_key, S.Zero) + S.One
        return dct

    @classmethod
    def transition_fundamental_slide(cls, key):
        from schubmult.abc import x
        from schubmult.rings.polynomial_algebra.fundamental_slide_poly_basis import FundamentalSlidePolyBasis

        key = tuple(key)
        length = len(key)
        total = sum(key)
        poly_basis = FundamentalSlidePolyBasis(x)

        ret = {}
        for candidate in cls._weak_compositions(length, total):
            coeff = poly_basis.to_monoms(candidate).get(key, S.Zero)
            if coeff != S.Zero:
                ret[candidate] = coeff
        return ret



    @classmethod
    def transition(cls, other_basis):
        from .forest_basis import ForestBasis
        from .fundamental_slide_basis import FundamentalSlideBasis
        from .j_basis import JBasis
        from .jt_basis import JTBasis
        from .key_basis import KeyBasis

        # from .key_basis import KeyBasis
        from .monomial_slide_basis import MonomialSlideBasis
        from .nelementary_basis import NElementaryBasis
        from .schubert_basis import SchubertBasis
        from .z_basis import ZBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        if other_basis == WordBasis:
            return lambda x: {x: S.One}
        if other_basis == KeyBasis:
            return lambda x: cls.transition_key(x)
        if other_basis == NElementaryBasis:
            return lambda x: cls.transition_nelementary(x)
        if other_basis == JBasis:
            return lambda x: cls.transition_jbasis(x)
        if other_basis == JTBasis:
            return lambda x: cls.transition_jtbasis(x)
        if other_basis == ZBasis:
            return lambda x: cls.transition_zbasis(x)
        if other_basis == FundamentalSlideBasis:
            return lambda x: cls.transition_fundamental_slide(x)
        if other_basis == ForestBasis:
            return lambda x: cls.transition_forest(x)
        if other_basis == MonomialSlideBasis:
            return lambda x: cls.transition_monomial_slide(x)
        # if other_basis == KeyBasis:
        #     return lambda x: cls.transition_key(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))
