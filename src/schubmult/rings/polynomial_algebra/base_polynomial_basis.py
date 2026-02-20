from abc import ABC, abstractmethod

from schubmult.symbolic import S
from schubmult.utils.perm_utils import add_perm_dict


class PolynomialBasis(ABC):
    @property
    @abstractmethod
    def genset(self):
        raise NotImplementedError

    @abstractmethod
    def is_key(self, x):
        raise NotImplementedError

    @abstractmethod
    def as_key(self, x):
        raise NotImplementedError

    def attach_key(self, dct):
        return {self.as_key(k): v for k, v in dct.items()}

    @property
    @abstractmethod
    def zero_monom(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def monomial_basis(self):
        raise NotImplementedError

    @abstractmethod
    def transition(self, other_basis):
        raise NotImplementedError

    @abstractmethod
    def printing_term(self, k):
        raise NotImplementedError

    @staticmethod
    def compose_transition(tkeyfunc, output):
        ret = {}
        for key, v in output.items():
            ret = add_perm_dict(ret, {k: v * v0 for k, v0 in tkeyfunc(key).items()})
        return ret

    def change_tensor_basis(self, tensor_elem, basis1, basis2):
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
        return self.monomial_basis.expand(self.transition(self.monomial_basis)(dct))

    def coproduct(self, key):
        from ...utils._mul_utils import _tensor_product_of_dicts_first

        monom_version = self.transition(self.monomial_basis)({key: S.One})
        ret_dict = {}
        for k, v in monom_version.items():
            this_coproduct = {k0: v0 * v for k0, v0 in self.monomial_basis.coproduct(k).items()}
            for (mk1, mk2), v2 in this_coproduct.items():
                ret_dict = add_perm_dict(
                    ret_dict,
                    _tensor_product_of_dicts_first(self.monomial_basis.transition({mk1: v2}), self.monomial_basis.transition({mk2: S.One})),
                )
        return ret_dict

    def product(self, key1, key2, coeff=S.One):
        from .monomial_basis import MonomialBasis

        mnb = MonomialBasis(genset=self.genset)
        left = self.transition(mnb)({key1: coeff})
        right = self.transition(mnb)({key2: S.One})

        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, mnb.transition(self)(mnb.product(key_schub_left, key_schub_right, v * v2)))
        return ret
