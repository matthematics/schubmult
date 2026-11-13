"""Quantum (single) Schubert polynomial ring: the ``QSx`` interface.

`QuantumSingleSchubertRing` is a `QuantumDoubleSchubertRing` with a zero
coefficient alphabet. This module also re-exports the quantum double and
parabolic quantum rings (``QDSx``, ``QPSx``, ``QPDSx``) for convenience.
"""

from functools import cache

import schubmult.rings.schubert.schubert_ring as spr
from schubmult.combinatorics.permutation import Permutation
from schubmult.symbolic import Add, Mul, Pow, S, sympify
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
    """The ring of quantum Schubert polynomials ``S^q_w(x)``; ``QSx`` is the standard instance."""

    def __init__(self, genset, q_var=q_var):
        """Use ``q_var`` as the quantum alphabet (default: ``GeneratingSet("q")``)."""
        super().__init__(genset, poly_genset(0), q_var=q_var)

    def __hash__(self):
        return hash((self.genset, self.coeff_genset, self.q_var, "QBS"))

    def quantize(self, poly):
        """Quantize a polynomial: expand in classical Schubert polynomials, reinterpret each ``S_w`` as
        the quantum ``S^q_w``, and expand back to a polynomial.
        """
        r = spr.SingleSchubertRing(self.genset)
        ssp = r.from_expr(poly)
        sspq = self.from_dict(dict(ssp.items()))
        return sspq.as_polynomial()

    def _coerce_mul(self, other):
        """Accept elements of this ring or a quantum double ring over the same ``x`` alphabet."""
        if other.ring == self:
            return other
        if type(other.ring) is QuantumDoubleSchubertRing:
            if self.genset == other.ring.genset and self.q_var == other.ring.q_var:
                return other
        return None

    @cache
    def cached_product(self, u, v, basis2):
        """Structure constants: ``schubmult_q_fast`` when ``basis2`` is this ring, else ``schubmult_q_double_fast``."""
        if self == basis2:
            return self.single_mul({u: S.One}, v)
        return self.double_mul({u: S.One}, v, self.coeff_genset, basis2.coeff_genset)

    @cache
    def cached_positive_product(self, u, v, basis2):
        """Same as ``cached_product``."""
        return self.cached_product(u, v, basis2)

    def mul_expr(self, elem, x):
        """Multiply by an expression: single ``x`` variables via ``mult_poly_q``, ``Add``/``Mul``/``Pow``
        recursively, anything else as a coefficient.
        """
        x = sympify(x)
        _Add = Add
        _Mul = Mul
        _Pow = Pow
        ind = self.genset.index(x)
        if ind != -1:
            return self.from_dict(self.mult_poly_single(elem, x, self.genset))
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
        """Build an element from a permutation/Lehmer list, a classical or parabolic element (converted
        to the quantum basis), or a polynomial expression.
        """
        genset = self.genset
        if not isinstance(genset, GeneratingSet_base):
            raise TypeError
        if isinstance(x, list) or isinstance(x, tuple):
            elem = self.from_dict({Permutation(x): self.domain.one})
        elif isinstance(x, Permutation):
            elem = self.from_dict({x: self.domain.one})
        elif isinstance(x, spr.DoubleSchubertElement):
            if x.ring.genset == self.genset:
                return self.from_expr(x.as_polynomial())
        elif isinstance(x, ParabolicQuantumDoubleSchubertElement):
            return self.from_expr(x.as_polynomial())
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
