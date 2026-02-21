from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class FundamentalSlideBasis(FreeAlgebraBasis):
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
    def printing_term(cls, k):
        return GenericPrintingTerm(k, "FS")


    @classmethod
    def transition_schubert(cls, key):
        from ..combinatorial import RCGraphRing
        r = RCGraphRing()
        dct = {}
        all_rcs = r.monomial(*key)
        for rc in all_rcs:
            if rc.is_quasi_yamanouchi:
                dct[(rc.perm, len(rc))] = dct.get((rc.perm, len(rc)), S.Zero) + S.One
        return dct


    @classmethod
    def transition(cls, other_basis):
        from .schubert_basis import SchubertBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))
