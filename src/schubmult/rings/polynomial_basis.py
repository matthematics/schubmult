from functools import cached_property
from itertools import zip_longest

# If you use symbolic coefficients like S.One, import S from your symbolic module
from schubmult.perm_lib import Permutation
from schubmult.rings.abstract_schub_poly import GenericPrintingTerm
from schubmult.symbolic import Add, Mul, S

# Utility for combining dictionaries with permutation keys
from schubmult.utils.perm_utils import add_perm_dict

# If you use _tensor_product_of_dicts_first in coproduct
# If you use self.monomial_basis and transition_schubert, import them from their respective modules
# If transition_schubert is a method of PolynomialBasis or another class, import accordingly
# For example, if it's a method of PolynomialBasis, you don't need to import it separately
# If you use TensorRing in change_tensor_basis
from .schubert_ring import Sx


class PolynomialBasis:
    @property
    def numvars(self): ...

    @property
    def genset(self): ...

    def is_key(self, x): ...

    def as_key(self, x): ...

    def with_numvars(self, numvars):
        return self.__class__(numvars=numvars, genset=self.genset)

    @property
    def zero_monom(self): ...

    @property
    def monomial_basis(self): ...

    # @classmethod
    # def product(self, key1, key2, coeff=S.One): ...

    # @classmethod
    # def coproduct(self, key, coeff=S.One): ...

    def transition(self, other_basis): ...

    def printing_term(self, k): ...

    # @property
    # def zero_monom(self): ...

    @staticmethod
    def compose_transition(tkeyfunc, output):
        ret = {}
        for key, v in output.items():
            ret = add_perm_dict(ret, {k: v * v0 for k, v0 in tkeyfunc(key).items()})
        return ret

    def change_tensor_basis(self, tensor_elem, basis1, basis2):
        from .tensor_ring import TensorRing

        ring1 = tensor_elem.ring.rings[0]
        ring2 = tensor_elem.ring.rings[1]
        Tring2 = TensorRing(ring1.__class__(basis=basis1), ring2.__class__(basis=basis2))
        res = Tring2.zero
        for (key1, key2), v in tensor_elem.items():
            new_elem1 = ring1(*key1).change_basis(basis1)
            new_elem2 = ring2(*key2).change_basis(basis2)
            res += v * Tring2.ext_multiply(new_elem1, new_elem2)
        return res

    # @classmethod
    # @cache
    # def coproduct(self, key):
    #     from ._mul_utils import _tensor_product_of_dicts_first

    #     return PolynomialBasis.compose_transition(
    #         lambda x: _tensor_product_of_dicts_first(self.monomial_basis.transition(self)(x[0]), self.monomial_basis.transition(self)(x[1])),
    #         PolynomialBasis.compose_transition(lambda y: self.monomial_basis.coproduct(y), self.transition(self.monomial_basis)(key)),
    #     )

    # @classmethod
    # @cache
    def coproduct(self, key):
        from ._mul_utils import _tensor_product_of_dicts_first

        monom_version = self.transition(self.monomial_basis)({key: S.One})
        ret_dict = {}
        for k, v in monom_version.items():
            this_coproduct = {k0: v0 * v for k0, v0 in self.monomial_basis.coproduct(k).items()}
            for (mk1, mk2), v2 in this_coproduct.items():
                ret_dict = add_perm_dict(
                    ret_dict,
                    _tensor_product_of_dicts_first(self.monomial_basis.transition(self.with_numvars(len(mk1)))({mk1: v2}), self.monomial_basis.transition(self.with_numvars(len(mk2)))({mk2: S.One})),
                )
        return ret_dict

    def product(self, key1, key2, coeff=S.One):
        mnb = MonomialBasis(self.numvars, self.genset)
        left = self.transition(mnb)(key1)
        right = self.transition(mnb)(key2)
        ret = {}

        for key_schub_right, v in right.items():
            for key_schub_left, v2 in left.items():
                ret = add_perm_dict(ret, PolynomialBasis.compose_transition(mnb.transition(self), MonomialBasis.product(key_schub_left, key_schub_right, v * v2 * coeff)))
        return ret


class MonomialBasis(PolynomialBasis):
    def is_key(self, x):
        return isinstance(x, tuple | list) and len(x) <= self.numvars

    def as_key(self, x):
        if len(x) > self.numvars:
            raise ValueError(f"Key {x} exceeds the number of variables {self.numvars}.")
        return (*x, *([0] * (self.numvars - len(x))))

    def printing_term(self, k):
        return GenericPrintingTerm(str(self.expand_monom(k)), "")

    def coproduct(self, key):
        result_dict = {}
        key = self.as_key(key)
        for i in range(self.numvars + 1):
            result_dict[(key[:i], key[i:])] = S.One
        return result_dict

    @property
    def monomial_basis(self):
        return self

    @property
    def numvars(self):
        return self._numvars

    @property
    def genset(self):
        return self._genset

    def __init__(self, genset, numvars):
        self._numvars = numvars
        self._genset = genset

    def product(self, key1, key2, coeff=S.One):
        return {tuple(a + b for a, b in zip_longest(key1, key2, fillvalue=0)): coeff}

    def expand_monom(self, monom):
        return Mul(*[self.genset[i + 1] ** monom[i] for i in range(len(monom))])

    def expand(self, dct):
        return Add(*[v * self.expand_monom(k) for k, v in dct.items()])

    def transition(self, other_basis):
        if isinstance(other_basis, MonomialBasis):
            return lambda x: x
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: {other_basis.as_key(k): v for k, v in other_basis.ring.from_expr(Add(*[v * self.expand_monom(k) for k, v in x.items()])).items()}
        return None

    def from_expr(self, expr):
        from .variables import genset_dict_from_expr

        return {self.as_key(k): v for k, v in genset_dict_from_expr(expr, self.genset).items()}

    @property
    def zero_monom(self):
        return self.as_key([])


class SchubertPolyBasis(PolynomialBasis):
    def __hash__(self):
        return hash(self.numvars, self.ring)

    def with_numvars(self, numvars):
        return self.__class__(numvars=numvars, ring=self.ring)

    @property
    def numvars(self):
        return self._numvars

    @property
    def genset(self):
        return self.ring.genset

    def printing_term(self, k):
        return self.ring.printing_term(k)

    def is_key(self, x):
        return isinstance(x, list | tuple | Permutation)

    def as_key(self, x):
        return (Permutation(x), self.numvars)

    def __init__(self, numvars, ring=None):
        self._numvars = numvars
        self.ring = ring
        if self.ring is None:
            self.ring = Sx([]).ring

    @cached_property
    def monomial_basis(self):
        return MonomialBasis(numvars=self.numvars, genset=self.ring.genset)

    def product(self, key1, key2, coeff=S.One):
        if key1[1] != key2[1] or key1[1] != self.numvars or key2[1] != self.numvars:
            return {}
        return {(k, self.numvars): v for k, v in self.ring.mul(self.ring.from_dict({key1[0]: coeff}), self.ring(key2[0])).items()}

    @property
    def zero_monom(self):
        return self.as_key([])

    def from_expr(self, expr):
        return self.ring.from_expr(expr)

    def transition(self, other_basis):
        if isinstance(other_basis, SchubertPolyBasis):
            return lambda x: x
        # if isinstance(other_basis, EXBasis):
        #     return lambda x: {S.One: self.ring(x).as_polynomial()}
        if isinstance(other_basis, MonomialBasis):
            from .variables import genset_dict_from_expr

            return lambda x: genset_dict_from_expr(self.ring.from_dict({k[0]: v for k, v in x.items()}).as_polynomial(), other_basis.genset)
        return None
