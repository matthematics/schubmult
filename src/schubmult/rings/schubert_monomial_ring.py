"""
Schubert Monomial Ring module

Provides base classes for rings whose basis elements represent Schubert monomials
(e.g., RC-graphs, BPDs, pipe dreams) with common operations like expansion to
polynomials, divided differences, and crystal operations.
"""

from __future__ import annotations

from schubmult.rings.abstract_schub_poly import TypedPrintingTerm
from schubmult.rings.base_schubert_ring import BaseSchubertElement, BaseSchubertRing
from schubmult.symbolic import S, sympify, sympy_Mul


class SchubertMonomialPrintingTerm(TypedPrintingTerm):
    """
    Printing term for Schubert monomial basis elements.

    Delegates printing to the underlying key object (typically an RCGraph, BPD, etc.)
    """



class SchubertMonomialRingElement(BaseSchubertElement):
    """
    Base class for ring elements whose basis elements are Schubert monomials.

    This provides a common interface for objects like:
    - RCGraphRingElement (basis elements are RCGraphs)
    - BPDRingElement (basis elements are BPDs)

    Common operations include:
    - Polynomial expansion via polyvalue()
    - Divided difference operators
    - Crystal structure operations (if the basis elements support them)
    """

    def polyvalue(self, x, y=None, **kwargs):
        """
        Evaluate as a polynomial in variables x (and optionally y).

        Linear extension: for each basis element, call its polyvalue() method
        and sum the results weighted by coefficients.

        Args:
            x: Variable or sequence of variables for polynomial evaluation
            y: Optional second set of variables for double Schubert polynomials
            **kwargs: Additional arguments passed to basis element polyvalue

        Returns:
            Symbolic expression representing the polynomial
        """
        from schubmult.symbolic import S

        result = S.Zero
        for basis_elem, coeff in self.items():
            result += coeff * basis_elem.polyvalue(x, y, **kwargs)
        return result

    def as_ordered_terms(self, *_, **__):
        if len(self.keys()) == 0:
            return [S.Zero]
        return [self[k] if k == self.ring.zero_monom else sympy_Mul(self[k], self.ring.printing_term(k)) for k in self.keys()]

    def to_free_algebra_element(self):
        """
        Convert to FreeAlgebra element in Schubert basis.
        """
        from schubmult.rings.free_algebra import FreeAlgebra
        from schubmult.rings.free_algebra_basis import SchubertBasis
        ASx = FreeAlgebra(SchubertBasis)
        ret = ASx.zero
        for monom, coeff in self.items():
            ret += coeff * ASx(monom.perm, len(monom))
        return ret


class SchubertMonomialRing(BaseSchubertRing):
    """
    Base class for rings whose basis elements are Schubert monomials.

    Inherits from BaseSchubertRing to provide standard ring operations (add, sub, mul, etc.)
    """

    def printing_term(self, key):
        return SchubertMonomialPrintingTerm(key)

    def from_dict(self, dct):
        elem = self.dtype()
        elem.update(dct)
        return elem

    def mul(self, a, b):
        if isinstance(b, SchubertMonomialRingElement):
            result_dict = {}
            for g1, c1 in a.items():
                for g2, c2 in b.items():
                    # RCGraph.product returns a dict {RCGraph: coeff}
                    prod = g1.product(g2)
                    for g3, c3 in prod.items():
                        result_dict[g3] = result_dict.get(g3, 0) + c1 * c2 * c3
            return self.from_dict(result_dict)
        try:
            result_dict = {k: v * sympify(b) for k, v in a.items()}
            return self.from_dict(result_dict)
        except Exception:
            raise NotImplementedError(f"Multiplication with {type(b)} not implemented for SchubertMonomialRingElement")

    @property
    def zero(self):
        return self.dtype()

    def __call__(self, key):
        return self.from_dict({key: 1})
