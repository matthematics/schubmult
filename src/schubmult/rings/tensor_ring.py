from functools import cache

from schubmult.symbolic import Mul, S, sympy_Mul
from schubmult.utils.logging import get_logger
from schubmult.utils.perm_utils import add_perm_dict

from ._mul_utils import _tensor_product_of_dicts, _tensor_product_of_dicts_first
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
        self._rings = rings
        genset = set()
        for r in self._rings:
            genset.update(set(r.genset))
        genset = tuple(genset)
        coeff_genset = set()
        for r in self._rings:
            if r.coeff_genset.label:
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
                dct = self.rings[0].from_dict({k1[0]: v1 * v2}) * elem2.ring.rings[0].from_dict({k2[0]: 1})
                for i in range(1, len(self.rings)):
                    dct2 = self.rings[i].from_dict({k1[i]: 1}) * elem2.ring.rings[i].from_dict({k2[i]: 1})
                    if i == 1:
                        dct = _tensor_product_of_dicts_first(dct, dct2)
                    else:
                        dct = _tensor_product_of_dicts(dct, dct2)
                ret_dict = add_perm_dict(ret_dict, dct)
        return self.from_dict(ret_dict)

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):
        # have to pull out gens
        if not isinstance(x.ring, TensorRing):
            if set(x.ring.genset) == set(self.genset):
                return x.coproduct(*[x.ring.genset.index(v) for v in self.rings[0].genset[1:]])
        return None

    @property
    def one(self):
        return self.from_dict({self.zero_monom: S.One})

    @cache
    def cached_schubpoly(self, k):
        return Mul(*[self.rings[i].cached_schubpoly(k[i]) for i in range(len(self.rings))])

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

    # def from_sympy(self, x):
    #     dct = {}
    #     elem1 = self.rings[0].from_sympy(x)
    #     dct = {(k,): v for k,v in elem1.items()}
    #     for i in range(1, len(self.rings)):
    #         dct_new = {}
    #         for k, v in dct.items():
    #             dct_new = add_perm_dict(dct_new,{(*k, k1): v1 for k1, v1 in self.rings[i].from_sympy(v).items()})
    #         dct = dct_new
    #     return self.from_dict(dct)

    def __call__(self, x):
        if isinstance(x, tuple):
            return self.from_dict({x: self.domain.one})
        return self.from_expr(x)


class TensorBasisElement(AbstractSchubPoly):
    is_commutative = False
    precedence = 50

    def __new__(cls, k, basis):
        return TensorBasisElement.__xnew_cached__(cls, k, basis)

    @staticmethod
    def __xnew__(_class, k, basis):
        obj =  AbstractSchubPoly.__new__(_class, k, None, None)
        obj._key = k
        obj.ring = basis
        return obj

    def __hash__(self):
        return hash((self._key, self.ring))

    @staticmethod
    @cache
    def __xnew_cached__(_class, k, basis):
        return TensorBasisElement.__xnew__(_class, k, basis)

    def _sympystr(self, printer):
        return " # ".join([printer._print(self.ring.rings[i].printing_term(self._key[i])) for i in range(len(self._key))])

    def _pretty(self, printer):
        return printer._print_TensorProduct(sympy_Mul(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))

    def _latex(self, printer):
        return printer._print_TensorProduct(sympy_Mul(*[self.ring.rings[i].printing_term(self._key[i]) for i in range(len(self._key))]))


class TensorRingElement(BaseSchubertElement):
    @property
    def free_symbols(self):
        ret = set()
        for k, v in self.items():
            ret.update(v.free_symbols)
            for i in range(len(k)):
                ret.update(self.ring.rings[i](k[i]).free_symbols)
        return ret

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [
            self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k))
            for k in self.keys()#sorted(self.keys(), key=lambda kkt: [(kk.inv, tuple(kk)) for kk in kkt])
        ]
