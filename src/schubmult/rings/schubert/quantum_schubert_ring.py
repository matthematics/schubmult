from functools import cache

import schubmult.mult.quantum as py
import schubmult.mult.quantum_double as yz
import schubmult.rings.schubert.schubert_ring as spr
from schubmult.combinatorial_reps.permutation import Permutation
from schubmult.symbolic import Add, Mul, Pow, S, sympify
from schubmult.symbolic.poly.poly_lib import xreplace_genvars
from schubmult.symbolic.poly.variables import GeneratingSet, GeneratingSet_base, poly_genset

from .parabolic_quantum_double_schubert_ring import (
    ParabolicQuantumDoubleSchubertElement,
    ParabolicQuantumDoubleSchubertRing,
    QPDSx,
    QPDSx_index,
    make_parabolic_quantum_basis,
)
from .parabolic_quantum_schubert_ring import ParabolicQuantumSingleSchubertRing, QPSx, make_single_parabolic_quantum_basis
from .quantum_double_schubert_ring import (
    QDSx,
    QuantumDoubleSchubertElement,
    QuantumDoubleSchubertPolynomial,
    QuantumDoubleSchubertRing,
    q_var,
)


class QuantumSingleSchubertRing(QuantumDoubleSchubertRing):
    def __init__(cls, genset):
        super().__init__(genset, poly_genset(0))

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, "QBS"))

    def quantize(self, poly):
        r = spr.SingleSchubertRing(self.genset)
        ssp = r.from_expr(poly)
        sspq = self.from_dict(dict(ssp.items()))
        return sspq.as_polynomial()

    def _coerce_mul(self, other):
        if other.ring == self:
            return other
        if type(other.ring) is QuantumDoubleSchubertRing:
            if self.genset == other.ring.genset:
                return other
        return None

    @cache
    def cached_product(self, u, v, basis2):
        if self == basis2:
            return py.schubmult_q_fast({u: S.One}, v)
        return {k: xreplace_genvars(x, self.coeff_genset, basis2.coeff_genset) for k, x in yz.schubmult_q_double_pair_generic(u, v).items()}

    @cache
    def cached_positive_product(self, u, v, basis2):
        return self.cached_product(u, v, basis2)

    def mul_expr(self, elem, x):
        x = sympify(x)
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(py.mult_poly_q(elem, x, self.genset))
        if isinstance(x, _Add):
            return self.sum([self.mul_expr(elem, arg) for arg in x.args])
        if isinstance(x, _Mul):
            res = elem
            for arg in x.args:
                res = self.mul_expr(res, arg)
            return res
        if isinstance(x, _Pow):
            res = elem
            for _ in range(int(x.args[1])):
                res = self.mul_expr(res, x.args[0])
            return res
        return self.from_dict({k: v * self.domain_new(x) for k, v in elem.items()})

    def new(self, x):
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, spr.DoubleSchubertElement):
            if x.genset == self.genset:
                return x.as_quantum()
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return x.as_quantum()
        else:
            elem = self.from_expr(x)
        return elem


QSx = QuantumSingleSchubertRing(GeneratingSet("x"))

__all__ = [
    "ParabolicQuantumDoubleSchubertElement",
    "ParabolicQuantumDoubleSchubertRing",
    "ParabolicQuantumSingleSchubertRing",
    "QDSx",
    "QPDSx",
    "QPDSx_index",
    "QPSx",
    "QSx",
    "QuantumDoubleSchubertElement",
    "QuantumDoubleSchubertPolynomial",
    "QuantumDoubleSchubertRing",
    "QuantumSingleSchubertRing",
    "make_parabolic_quantum_basis",
    "make_single_parabolic_quantum_basis",
    "q_var",
]
