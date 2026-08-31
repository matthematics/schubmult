"""
Multiplication algorithms for Schubert polynomials.

This module provides kernels for computing products of Schubert polynomials
in various settings:

- schubmult_py: Ordinary (single) Schubert polynomial multiplication
- schubmult_double: Double Schubert polynomial multiplication
- schubmult_q: Quantum Schubert polynomial multiplication
- schubmult_q_double: Quantum double Schubert polynomial multiplication
- grothmult_double: Double Grothendieck multiplication by a degree-one class

Also includes positivity utilities (posify, compute_positive_rep) for root-based representations.
"""

from schubmult.mult.double import mult_poly_double, schubmult_double
from schubmult.mult.groth_double import (
    epsilon_chain,
    grothmult_double,
    monk_chain,
    mult_poly_groth_double,
    single_variable_groth,
)
from schubmult.mult.positivity import compute_positive_rep, posify
from schubmult.mult.quantum import schubmult_q
from schubmult.mult.quantum_double import factor_out_q, schubmult_q_double
from schubmult.mult.separated_descents import (
    grothmult_double as separated_descents_grothmult_double,
)
from schubmult.mult.separated_descents import separated_descents_coeffs
from schubmult.mult.single import mult_poly_py, schubmult_py

__all__ = [
    "compute_positive_rep",
    "epsilon_chain",
    "factor_out_q",
    "grothmult_double",
    "monk_chain",
    "mult_poly_double",
    "mult_poly_groth_double",
    "mult_poly_py",
    "posify",
    "schubmult_double",
    "schubmult_py",
    "schubmult_q",
    "schubmult_q_double",
    "separated_descents_coeffs",
    "separated_descents_grothmult_double",
    "single_variable_groth",
]
