from schubmult.combinatorics.indexed_forests import word_to_indexed_forest
from schubmult.rings.free_algebra.free_algebra_basis import FreeAlgebraBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S


class ForestBasis(FreeAlgebraBasis):
    @classmethod
    def is_key(cls, x):
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        return tuple(x)

    # @classmethod
    # def from_rc_graph(cls, rc_graph):
    #     return {rc_graph.length_vector(): 1}

    # @classmethod
    # def product(cls, Forest1, Forest2, coeff=S.One):
    #     return {(*Forest1, *Forest2): coeff}

    zero_monom = ()


    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "Forest")


    @classmethod
    def transition_schubert(cls, key):
        from schubmult.rings.combinatorial import RCGraphRing
        r = RCGraphRing()
        dct = {}
        all_rcs = r.monomial(*key)

        for rc in all_rcs:
            # if rc.is_lowest_weight:
            #     dct[(rc.perm, len(rc))] = dct.get((rc.perm, len(rc)), S.Zero) + S.One
            word = list(reversed(rc.perm_word))
            indfor = word_to_indexed_forest(word)
            if indfor.code == key:
                dct[(rc.perm, len(rc))] = dct.get((rc.perm, len(rc)), S.Zero) + S.One
        return dct


    @classmethod
    def transition(cls, other_basis):
        from schubmult.rings.free_algebra.schubert_basis import SchubertBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))

if __name__ == "__main__":
    from schubmult import FreeAlgebra, SchubertBasis

    Forest = FreeAlgebra(ForestBasis)
    print(Forest((0, 2, 0, 1)).change_basis(SchubertBasis))
