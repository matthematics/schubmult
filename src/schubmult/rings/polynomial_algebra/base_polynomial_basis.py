from abc import ABC, abstractmethod

from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict


class PolynomialBasis(ABC):
    """Abstract base class for polynomial algebra bases.

    Subclasses define how keys are represented, how to transition between
    bases, and how to expand elements into explicit polynomials. Default
    implementations delegate through the :class:`MonomialBasis`.
    """
    @property
    def genset(self):
        return self._genset

    @abstractmethod
    def is_key(self, x):
        """Return True if *x* is a valid key for this basis."""
        raise NotImplementedError

    @abstractmethod
    def as_key(self, x):
        """Normalize *x* into a canonical key for this basis."""
        raise NotImplementedError

    def attach_key(self, dct):
        """Normalize all keys in *dct* via :meth:`as_key`."""
        return {self.as_key(k): v for k, v in dct.items()}

    @property
    @abstractmethod
    def zero_monom(self):
        raise NotImplementedError

    @property
    def monomial_basis(self):
        return self._monomial_basis

    @abstractmethod
    def transition(self, other_basis):
        """Return a function mapping dicts of this basis to dicts in *other_basis*."""
        raise NotImplementedError

    def from_expr(self, expr, length=None):
        """Parse a symbolic expression into this basis."""
        from .monomial_basis import MonomialBasis

        monomial_basis = MonomialBasis(genset=self.genset)
        monomial_dict = monomial_basis.from_expr(expr, length=length)
        if isinstance(self, MonomialBasis):
            return monomial_dict
        return monomial_basis.transition(self)(monomial_dict)

    @abstractmethod
    def printing_term(self, k):
        """Return the display symbol for key *k*."""
        raise NotImplementedError

    @staticmethod
    def compose_transition(tkeyfunc, output):
        """Apply a transition function to a dict of basis elements."""
        return tkeyfunc(output)

    @classmethod
    def change_tensor_basis(cls, tensor_elem, basis1, basis2):
        """Change the bases of both factors of a tensor element.

        Args:
            tensor_elem: An element of a tensor product ring.
            basis1: Target basis for the left factor.
            basis2: Target basis for the right factor.

        Returns:
            The tensor element re-expressed in the new bases.
        """
        from ..tensor_ring import TensorRing

        ring1 = tensor_elem.ring.rings[0]
        ring2 = tensor_elem.ring.rings[1]
        Tring2 = TensorRing(ring1.__class__(basis=basis1), ring2.__class__(basis=basis2))
        res = Tring2.zero
        for (key1, key2), v in tensor_elem.items():
            new_elem1 = ring1(*key1).change_basis(basis1)
            new_elem2 = ring2(*key2).change_basis(basis2)
            res += v * Tring2.ext_multiply(new_elem1, new_elem2)
        return res

    def expand(self, dct):
        """Expand a basis dict into an explicit polynomial expression."""
        return self.monomial_basis.expand(self.transition(self.monomial_basis)(dct))

    def coproduct(self, key):
        """Compute the coproduct of *key* by delegating through the monomial basis."""
        from ...utils._mul_utils import _tensor_product_of_dicts_first

        monom_version = self.transition(self.monomial_basis)({key: S.One})
        ret_dict = {}
        for k, v in monom_version.items():
            this_coproduct = {k0: v0 * v for k0, v0 in self.monomial_basis.coproduct(k).items()}
            for (mk1, mk2), v2 in this_coproduct.items():
                ret_dict = add_perm_dict(
                    ret_dict,
                    _tensor_product_of_dicts_first(self.monomial_basis.transition(self)({mk1: v2}), self.monomial_basis.transition(self)({mk2: S.One})),
                )
        return ret_dict

    def product(self, key1, key2, coeff=S.One):
        """Multiply two keys by transitioning to the monomial basis and back."""
        from .monomial_basis import MonomialBasis

        mnb = MonomialBasis(genset=self.genset)
        left = self.transition(mnb)({key1: coeff})
        right = self.transition(mnb)({key2: S.One})

        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, mnb.transition(self)(mnb.product(key_schub_left, key_schub_right, v * v2)))
        return ret

    def __init__(self, genset, *args, **kwargs):  # noqa: ARG002
        self._genset = genset

    @classmethod
    def dual_basis(cls):
        """Return the dual free algebra basis class."""
