from functools import cache

import numpy as np
from sympy import Tuple

from schubmult.symbolic import Mul, S, sympy_Mul
from schubmult.utils.logging import get_logger

from .abstract_schub_poly import AbstractSchubPoly
from .base_schubert_ring import BaseSchubertElement, BaseSchubertRing

logger = get_logger(__name__)


class ProductRing(BaseSchubertRing):
    # tensor ring
    def __eq__(self, other):
        return type(self) is type(other) and self.rings == other.rings

    @property
    def args(self):
        return ()



    # def from_rc_graph_tensor(self, rc_graph_tensor):
    #     return self.ext_multiply(self.rings[0].from_rc_graph(rc_graph_tensor[0]), self.rings[1].from_rc_graph(rc_graph_tensor[1]))

    # def sub(self, elem, other):
    #     print("Mooger")
    #     return self.from_dict(add_perm_dict(elem, {k: -v for k,v in other.items()}))

    def __init__(self, *rings):
        self._rings = rings
        ring_list = [*self._rings]
        while any(isinstance(r, ProductRing) for r in ring_list):
            new_ring_list = []
            for r in ring_list:
                if isinstance(r, ProductRing):
                    new_ring_list.extend(r.rings)
                else:
                    new_ring_list.append(r)
            ring_list = new_ring_list
        self._rings = tuple(ring_list)
        genset = set()
        for r in self._rings:
            try:
                genset.update(set(r.genset))
            except AttributeError:
                pass
        genset = tuple(genset)
        coeff_genset = set()
        for r in self._rings:
            try:
                if r.coeff_genset.label:
                    coeff_genset.update(set(r.coeff_genset))
            except AttributeError:
                pass
        coeff_genset = tuple(coeff_genset)
        super().__init__(list(genset), list(coeff_genset))
        self.zero_monom = tuple([self.rings[i].zero_monom for i in range(len(self.rings))])
        # self.dtype = type("ProductRingElement", (ProductRingElement,), {"ring": self})

    def dtype(self):
        elem = ProductRingElement()
        elem.ring = self
        return elem

    def __hash__(self):
        return hash(self.rings)

    @property
    def rings(self):
        return self._rings

    def rmul(self, elem1, elem2):
        # print(f"{dict(elem1)=} {elem2=}")
        # print(f"{self.zero_monom=} {type(elem2)=}")
        # return self.from_dict({self.zero_monom: elem2}) * elem1
        return self.from_dict({k: v * elem2 for k, v in elem1.items()})
        # except Exception:
        #     # import traceback
        #     # traceback.print_exc()
        #     raise

    def from_dict(self, element):
        dct = self.dtype()
        dct.update({k: v for k, v in element.items() if v != 0})
        return dct

    def mul(self, elem1, elem2):
        return self.from_dict({k: v * elem2 for k, v in elem1.items()})

    def _coerce_add(self, x):  # noqa: ARG002
        return None

    def _coerce_mul(self, x):
        # have to pull out gens
        if not isinstance(x.ring, ProductRing):
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
        return ProductBasisElement(k, self)


    def new(self, x):
        elem = self.dtype()
        elem._arr = np.array(x, dtype=object)
        return elem

    def __call__(self, *x):
        if len(x) == 1 and isinstance(x[0], self.dtype()):
            return self.new(x[0]._arr.copy())
        if len(x) == 1 and isinstance(x[0], list | tuple | np.ndarray):
            return self.new([self.rings[i](x[0][i]) for i in range(len(self.rings))])
        return self.new([self.rings[i](x[i]) for i in range(len(self.rings))])



class ProductBasisElement(AbstractSchubPoly):
    is_commutative = False
    precedence = 50

    def __init__(self, elem):
        self._arr = np.array(elem._arr, dtype=object)

    def __hash__(self):
        return hash((self._arr.tobytes(), self.ring))

    def _sympystr(self, printer):
        return " # ".join([printer._print(elem.ring.printing_term(elem)) for elem in self._arr])

    def _pretty(self, printer):
        return printer._print_seq(Tuple(*[self.ring.rings[i].printing_term(self._arr[i]) for i in range(len(self._arr))]))

    def _latex(self, printer):
        return printer._print_seq(Tuple(*[self.ring.rings[i].printing_term(self._arr[i]) for i in range(len(self._arr))]))


class ProductRingElement(BaseSchubertElement):
    def __init__(self):
        pass


    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [
            self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k))
            for k in self.keys()  # sorted(self.keys(), key=lambda kkt: [(kk.inv, tuple(kk)) for kk in kkt])
        ]

    def __add__(self, other):
        return self.ring(self._arr + other._arr)

    def __radd__(self, other):
        return self.ring(other + self._arr)

    def __sub__(self, other):
        return self.ring(self._arr - other._arr)

    def __rsub__(self, other):
        return self.ring(other - self._arr)

    def __neg__(self):
        return self.ring(np.array([-v for v in self._arr], dtype=object))

    def __pow__(self, val):
        if val == 0:
            return self.ring.one
        if val == 1:
            return self
        return self.ring(self._arr ** val)

    def __mul__(self, other):
        return self.ring(self._arr * other._arr)

    def __rmul__(self, other):
        return self.ring(other * self._arr)
