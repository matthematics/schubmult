from functools import cache

import sympy
from symengine import S

import schubmult.rings._utils as utils
from schubmult.perm_lib import Permutation
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from .abstract_schub_poly import AbstractSchubPoly
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing

logger = get_logger(__name__)


class TensorRing(BaseSchubertRing):
    # tensor ring
    def __eq__(self, other):
        return type(self) is type(other) and self.rings == other.rings

    @property
    def args(self):
        return ()

    def __init__(self, *rings):
        self._rings = list(rings)
        new_rings = []
        unpacked = True
        while unpacked:
            unpacked = False
            for r in self._rings:
                if isinstance(r, TensorRing):
                    unpacked = True
                    for r2 in r.rings:
                        if r2 not in new_rings:
                            new_rings += [r2]
                else:
                    if r not in new_rings:
                        new_rings += [r]
            self._rings = new_rings
        self._rings = tuple(new_rings)
        genset = set()
        for r in self._rings:
            genset.update(set(r.genset))
        genset = tuple(genset)
        coeff_genset = set()
        for r in self._rings:
            coeff_genset.update(set(r.coeff_genset))
        coeff_genset = tuple(coeff_genset)
        super().__init__(list(genset), list(coeff_genset))
        self.zero_monom = tuple([self.rings[i].zero_monom for i in range(len(self.rings))])
        self.dtype = type("TensorRingElement", (TensorRingElement,), {"ring": self})

    def __hash__(self):
        return hash(self._rings)

    @property
    def rings(self):
        return self._rings

    def mul(self, elem1, elem2):
        ret_dict = {}
        for k1, v1 in elem1.items():
            for k2, v2 in elem2.items():
                dct = self.rings[0].from_dict({k1[0]: v1 * v2}) * self.rings[0].from_dict({k2[0]: 1})
                for i in range(1, len(self.rings)):
                    dct2 = self.rings[i].from_dict({k1[i]: 1}) * self.rings[i].from_dict({k2[i]: 1})
                    dct = utils._tensor_product_of_dicts(dct, dct2)
                ret_dict = add_perm_dict(ret_dict, dct)
        return self.from_dict(ret_dict)

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):  # noqa: ARG002
        return None

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    @cache
    def cached_schubpoly(self, k):
        return sympy.Mul(*[self.rings[i].cached_schubpoly(k[i]) for i in range(len(self.rings))])

    def printing_term(self, k):
        return TensorBasisElement(k, self)

    def from_comp_ring(self, t):
        dct = {}
        for k, v in t.items():
            new_k = list(self.zero_monom)
            if isinstance(t.ring, TensorRing):
                for i in range(len(t.ring.rings)):
                    new_k[self.rings.index(t.ring.rings[i])] = k[i]
            else:
                new_k[self.rings.index(t.ring)] = k
            dct[tuple(new_k)] = v
        return self.from_dict(dct)

    def from_sympy(self, x):
        elem1 = self.rings[0].from_sympy(x)
        for i in range(1, len(self.rings)):
            res = self.zero
            for k, v in elem1.items():
                res += self.from_dict(utils._tensor_product_of_dicts({k: S.One}, self.rings[i].from_sympy(v)))
            elem1 = res
        return elem1

    def __call__(self, x):
        if isinstance(x, tuple):
            return self.from_dict({x: 1})
        return self.from_sympy(x)


class TensorBasisElement(AbstractSchubPoly):
    is_commutative = False

    def __new__(cls, k, basis):
        return TensorBasisElement.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj = AbstractSchubPoly.__new__(_class, k, basis)
        if not basis.is_implicit:
            obj.precedence = 50
        else:
            obj.precedence = 1000
        return obj

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return TensorBasisElement.__xnew__(_class, k, basis)

    def _sympystr(self, printer):
        if self.ring.is_implicit:
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
        return self.ring.tensor_symbol.join([printer._print(self.ring.rings[i].printing_term(self._key[i])) for i in range(len(self._key))])

    def _pretty(self, printer):
        if self.ring.is_implicit:
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
        return printer._print_TensorProduct(sympy.Mul(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))

    def _latex(self, printer):
        if self.ring.is_implicit:
            if self._key == self.ring.zero_monom:
                return printer._print(S.One)
            for i in range(len(self._key)):
                if self._key[i] != Permutation([]):
                    return printer._print(self.ring.rings[i].printing_term(self._key[i]))
        return printer._print_TensorProduct(sympy.Mul(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))


class TensorRingElement(BaseSchubertElement):
    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            for i in range(len(k)):
                ret.update(self.rings.rings[i](k[i]).free_symbols)
        return ret

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [sympy.sympify(S.Zero)]
        return [
            self.ring.domain.to_sympy(self[k]) if k == self.ring.zero_monom else sympy.Mul(self.ring.domain.to_sympy(self[k]), self.ring.printing_term(k))
            for k in sorted(self.keys(), key=lambda kkt: [(kk.inv, tuple(kk)) for kk in kkt])
        ]
