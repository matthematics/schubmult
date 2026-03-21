from schubmult.combinatorics.permutation import Permutation, uncode
from schubmult.symbolic import S
from schubmult.utils.tuple_utils import pad_tuple

from ..printing import GenericPrintingTerm
from .free_algebra_basis import FreeAlgebraBasis
from .schubert_basis import SchubertBasis


class CompositionSchubertBasis(FreeAlgebraBasis):
    """Schubert basis indexed by padded trimcode compositions.

    A key is a composition ``c`` and corresponds to Schubert key
    ``(uncode(c), len(c))``.
    """

    @classmethod
    def is_key(cls, x):
        """Return True if *x* is a tuple or list (composition)."""
        return isinstance(x, (tuple, list))

    @classmethod
    def as_schubert_key(cls, key):
        """Convert a composition key to a Schubert key ``(Permutation, length)``."""
        return (uncode(key), len(key))

    @classmethod
    def as_key(cls, key):
        """Normalize a key to a composition tuple.

        Accepts either a Schubert key ``(Permutation, int)`` or a raw tuple.
        """
        if isinstance(key, tuple) and len(key) == 2 and isinstance(key[1], int) and isinstance(key[0], (Permutation, list, tuple)):
            perm = Permutation(key[0])
            return pad_tuple(perm.trimcode, key[1])
        return tuple(key)

    @classmethod
    def from_rc_graph(cls, rc_graph):
        """Return the composition key for the given RC graph."""
        out = {}
        for key, coeff in SchubertBasis.from_rc_graph(rc_graph).items():
            out[cls.as_key(key)] = coeff
        return out

    @classmethod
    def inject(cls, key1, i, key2, coeff=S.One):
        """Inject *key2* into *key1* at position *i* using Schubert multiplication."""
        from schubmult.utils.perm_utils import mu_A

        from ..schubert.schubert_ring import Sx
        w0_12_cd = list(range(len(uncode(key1)) + len(uncode(key2)) - 1, 0, -1))
        w0_12 = uncode(w0_12_cd)
        the_range = list(range(len(w0_12_cd)))
        mu1 = mu_A(w0_12_cd, the_range[:i] + the_range[i + len(key2):])
        mu2 = mu_A(w0_12_cd, the_range[i:i + len(key2)])
        mul_perm1 = uncode(key1) * uncode(mu1)
        mul_perm2 = uncode(key2) * uncode(mu2)

        bigprod = Sx(mul_perm1) * Sx(mul_perm2)
        ret = {}
        target_inv = sum(key1) + sum(key2)
        for up_w, coeff0 in bigprod.items():
            w = up_w * (~w0_12)
            if w.inv == target_inv and len(w.trimcode) <= len(key1) + len(key2):
                new_key = cls.as_key((w, len(key1) + len(key2)))
                ret[new_key] = ret.get(new_key, S.Zero) + coeff * coeff0
        return ret

    @classmethod
    def product(cls, key1, key2, coeff=S.One):
        """Multiply two composition keys by delegating to SchubertBasis."""
        return {cls.as_key(k): v for k, v in SchubertBasis.product(cls.as_schubert_key(key1), cls.as_schubert_key(key2), coeff).items()}

    zero_monom = ()

    @classmethod
    def coproduct(cls, key):
        """Compute the coproduct by delegating to SchubertBasis."""
        return {cls.as_key(k): v for k, v in SchubertBasis.coproduct(cls.as_schubert_key(key)).items()}

    @classmethod
    def bcoproduct(cls, key):
        """Compute the bar-coproduct by delegating to SchubertBasis."""
        return {cls.as_key(k): v for k, v in SchubertBasis.bcoproduct(cls.as_schubert_key(key)).items()}

    @classmethod
    def internal_product(cls, key1, key2, coeff=S.One):
        """Compute the internal product by delegating to SchubertBasis."""
        return {cls.as_key(k): v for k, v in SchubertBasis.internal_product(cls.as_schubert_key(key1), cls.as_schubert_key(key2), coeff).items()}

    @classmethod
    def skew_element(cls, w, u, n):
        """Compute the skew element by delegating to SchubertBasis."""
        return {cls.as_key(k): v for k, v in SchubertBasis.skew_element(w, u, n).items()}


    @classmethod
    def dual_basis(cls):
        """Return the dual basis (delegates to SchubertBasis)."""
        return SchubertBasis.dual_basis()

    @classmethod
    def transition(cls, other_basis):
        """Return a transition function from CompositionSchubertBasis to *other_basis*."""
        if other_basis == CompositionSchubertBasis:
            return lambda x: {x: S.One}
        if other_basis == SchubertBasis:
            return lambda x: {cls.as_schubert_key(x): S.One}
        return lambda x: SchubertBasis.transition(other_basis)(cls.as_schubert_key(x))

    @classmethod
    def printing_term(cls, k):
        """Return a ``CompSchub``-labelled display object for the composition key *k*."""
        return GenericPrintingTerm(cls.as_key(k), "CompSchub")
