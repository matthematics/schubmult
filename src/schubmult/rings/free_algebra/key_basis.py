from schubmult.rings.free_algebra.free_algebra_basis import FreeAlgebraBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S


class KeyBasis(FreeAlgebraBasis):
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
    # def product(cls, key1, key2, coeff=S.One):
    #     return {(*key1, *key2): coeff}

    zero_monom = ()

    @classmethod
    def dual_basis(cls):
        from ..polynomial_algebra.key_poly_basis import KeyPolyBasis
        return KeyPolyBasis

    @classmethod
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "Key")


    @classmethod
    def transition_schubert(cls, key):
        from schubmult.rings.combinatorial import RCGraphRing
        r = RCGraphRing()
        dct = {}
        all_rcs = r.monomial(*key)
        for rc in all_rcs:
            if rc.extremal_weight == key:
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

    Key = FreeAlgebra(KeyBasis)
    print(Key((0, 2, 0, 1)).change_basis(SchubertBasis))
