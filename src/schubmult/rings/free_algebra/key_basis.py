from functools import cache

from schubmult.rings.free_algebra.free_algebra_basis import FreeAlgebraBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S


class KeyBasis(FreeAlgebraBasis):
    """Key polynomial (Demazure character) basis of the free algebra.

    Keys are weak composition tuples. Transitions to the Schubert basis
    use RC graph enumeration filtered by extremal weight.
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a tuple or list."""
        return isinstance(x, tuple | list)

    @classmethod
    def as_key(cls, x):
        """Normalize *x* to a tuple key."""
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
        """Return the KeyPolyBasis as the dual of KeyBasis."""
        from ..polynomial_algebra.key_poly_basis import KeyPolyBasis
        return KeyPolyBasis

    @classmethod
    def printing_term(cls, k):
        """Return a ``Key``-labelled display object for key *k*."""
        return GenericPrintingTerm(k, "Key")


    @classmethod
    def transition_schubert(cls, key):
        """Transition a key composition to the Schubert basis via RC graphs."""
        from schubmult import RCGraphRing
        #from schubmult import uncode
        r = RCGraphRing()
        dct = {}
        #all_rcs = {rc: coeff for rc, coeff in r.monomial(*key).items() if rc.extremal_weight == key}
        all_rcs = r.monomial(*key)
        #all_perms = {rc.perm for rc in all_rcs.keys()}
        #for perm in all_perms:
        for rc in all_rcs.keys():
            #for hw_rc in RCGraph.all_hw_rcs(perm, len(key)):
            if rc.extremal_weight == key:
                dct[(rc.perm, len(key))] = dct.get((rc.perm, len(key)), S.Zero) + 1
        # seen = set()
        # for rc, coeff in all_rcs.items():
        #     if rc.extremal_weight == key and rc.to_highest_weight()[0] not in seen:
        #         dct[(rc.perm, len(rc))] = dct.get((rc.perm, len(rc)), S.Zero) + S.One
        #         seen.add(rc.to_highest_weight()[0])
        return dct

    @classmethod
    @cache
    def product(cls, key1, key2, coeff=S.One):
        """Multiply two keys by transitioning to WordBasis and back."""
        from schubmult.utils._mul_utils import add_perm_dict

        from .schubert_basis import SchubertBasis

        left = cls.transition_schubert(key1)
        right = cls.transition_schubert(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, FreeAlgebraBasis.compose_transition(SchubertBasis.transition(cls), SchubertBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret


    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from KeyBasis to *other_basis*."""
        from schubmult.rings.free_algebra.schubert_basis import SchubertBasis

        if other_basis == SchubertBasis:
            return lambda x: cls.transition_schubert(x)
        return lambda x: FreeAlgebraBasis.compose_transition(SchubertBasis.transition(other_basis), cls.transition_schubert(x))

if __name__ == "__main__":
    from schubmult import FreeAlgebra, SchubertBasis

    Key = FreeAlgebra(KeyBasis)
    print(Key((0, 2, 0, 1)).change_basis(SchubertBasis))
