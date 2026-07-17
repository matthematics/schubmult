
from schubmult.symbolic import S

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis


class GlideBasis(FreeAlgebraBasis):
    """Glide basis of the free algebra.

    Keys are weak composition tuples. Transitions are computed via
    polynomial algebra slide polynomial expansions.
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
    def printing_term(cls, k):
        """Return an ``FS``-labelled display object for key *k*."""
        return GenericPrintingTerm(k, "FS")


    @classmethod
    def transition_grothendieck(cls, key):
        """Transition a glide key to the Grothendieck basis."""
        from schubmult.rings.combinatorial.wc_graph_ring import WCGraphRing

        wr = WCGraphRing()

        ring_result = wr.monomial(*key)
        dct = {}
        for wc, coeff in ring_result.items():
            if wc.is_quasi_yamanouchi:
                dct[(wc.perm, len(wc))] = dct.get((wc.perm, len(wc)), S.Zero) + coeff
        return dct


    @classmethod
    def dual_basis(cls):
        """Return the GlidePolyBasis as the dual."""
        from ..polynomial_algebra.glide_poly_basis import GlidePolyBasis
        return GlidePolyBasis

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from GlideBasis to *other_basis*."""
        from .grothendieck_basis import GrothendieckBasis

        if other_basis == GrothendieckBasis:
            return lambda x: cls.transition_grothendieck(x)
        return lambda x: FreeAlgebraBasis.compose_transition(GrothendieckBasis.transition(other_basis), cls.transition_grothendieck(x))
