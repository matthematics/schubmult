from schubmult.rings.free_algebra._core import FreeAlgebra
from schubmult.rings.free_algebra.free_algebra_basis import FreeAlgebraBasis
from schubmult.rings.printing import GenericPrintingTerm
from schubmult.symbolic import S


class GroveBasis(FreeAlgebraBasis):
    """Grove basis of the free algebra.

    Keys are tuples representing indexed-grove weight vectors.
    Transitions to the Grothendieck basis use RC graph enumeration
    filtered by grove weight.
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
    # def product(cls, Grove1, Grove2, coeff=S.One):
    #     return {(*Grove1, *Grove2): coeff}

    zero_monom = ()

    @classmethod
    def printing_term(cls, k):
        """Return a ``GroveDual``-labelled display object for key *k*."""
        return GenericPrintingTerm(k, "GroveDual")

    @classmethod
    def dual_basis(cls):
        """Return the GrovePolyBasis as the dual of GroveBasis."""
        from ..polynomial_algebra.grove_poly_basis import GrovePolyBasis

        return GrovePolyBasis

    @classmethod
    def transition_grothendieck(cls, key):
        """Transition a grove key to the Grothendieck basis via WC graph enumeration."""
        from schubmult.rings.combinatorial import WCGraphRing

        r = WCGraphRing()
        dct = {}
        all_rcs = r.monomial(*key)

        for wc in all_rcs:
            if wc.forest_weight == key:
                dct[(wc.perm, len(wc))] = dct.get((wc.perm, len(wc)), S.Zero) + S.One
        return dct

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from GroveBasis to *other_basis*."""
        from schubmult.rings.free_algebra.grothendieck_basis import GrothendieckBasis

        if other_basis == GrothendieckBasis:
            return lambda x: cls.transition_grothendieck(x)
        return lambda x: FreeAlgebraBasis.compose_transition(GrothendieckBasis.transition(other_basis), cls.transition_grothendieck(x))

GroveDual = FreeAlgebra(GroveBasis)


if __name__ == "__main__":
    from schubmult import FreeAlgebra, GrothendieckBasis

    Grove = FreeAlgebra(GroveBasis)
    print(Grove((0, 2, 0, 1)).change_basis(GrothendieckBasis))
