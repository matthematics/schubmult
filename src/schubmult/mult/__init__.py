"""
Multiplication algorithms for Schubert polynomials.

This module provides kernels for computing products of Schubert polynomials
in various settings:

- schubmult_py: Ordinary (single) Schubert polynomial multiplication
- schubmult_double: Double Schubert polynomial multiplication
- schubmult_q: Quantum Schubert polynomial multiplication
- schubmult_q_double: Quantum double Schubert polynomial multiplication

Also includes positivity display utilities for root-based representations.
"""

from schubmult.mult.double import mult_poly_double, schubmult_double
from schubmult.mult.positivity import display_positive, factor_out_q
from schubmult.mult.quantum import schubmult_q
from schubmult.mult.quantum_double import schubmult_q_double
from schubmult.mult.single import mult_poly_py, schubmult_py

__all__ = [
    "display_positive",
    "factor_out_q",
    "mult_poly_double",
    "mult_poly_py",
    "schubmult_double",
    "schubmult_py",
    "schubmult_q",
    "schubmult_q_double",
]
