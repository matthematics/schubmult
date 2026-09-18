"""
Multiplication algorithms for Schubert polynomials.

This module provides kernels for computing products of Schubert polynomials
in various settings:

- schubmult_py: Ordinary (single) Schubert polynomial multiplication
- schubmult_double: Double Schubert polynomial multiplication
- schubmult_q: Quantum Schubert polynomial multiplication
- schubmult_q_double: Quantum double Schubert polynomial multiplication
- grothmult_double: Double Grothendieck multiplication by a degree-one class
- grothmult_q_double: Quantum double Grothendieck multiplication (conjectural Molev--Sagan rule)
- grothmult_q: Quantum (single) Grothendieck multiplication

Also includes positivity utilities (posify, compute_positive_rep) for root-based representations.
"""

from schubmult.mult.double import mult_poly_double, schubmult_double
from schubmult.mult.groth import grothmult_py, mult_poly_groth
from schubmult.mult.groth_double import (
    elem_sym_perms_groth,
    epsilon_chain,
    groth_elem_sym_poly,
    grothmult_double,
    grothmult_double_block,
    grothmult_double_pieri,
    monk_chain,
    mult_poly_groth_double,
    one_plus_beta_x_groth,
    single_variable_groth,
)
from schubmult.mult.groth_quantum import grothmult_q, grothmult_q_dict, grothmult_q_pieri
from schubmult.mult.groth_quantum_double import (
    groth_elem_sym_poly_q,
    grothmult_q_double,
    grothmult_q_double_dict,
    grothmult_q_double_pieri,
    grothmult_q_double_top,
    lm_quantize,
    qgroth_poly,
    quantum_elem_sym,
    quantum_pieri_chains,
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
    "elem_sym_perms_groth",
    "epsilon_chain",
    "factor_out_q",
    "groth_elem_sym_poly",
    "groth_elem_sym_poly_q",
    "grothmult_double",
    "grothmult_double_block",
    "grothmult_double_pieri",
    "grothmult_py",
    "grothmult_q",
    "grothmult_q_dict",
    "grothmult_q_double",
    "grothmult_q_double_dict",
    "grothmult_q_double_pieri",
    "grothmult_q_double_top",
    "grothmult_q_pieri",
    "lm_quantize",
    "monk_chain",
    "mult_poly_double",
    "mult_poly_groth",
    "mult_poly_groth_double",
    "mult_poly_py",
    "one_plus_beta_x_groth",
    "posify",
    "qgroth_poly",
    "quantum_elem_sym",
    "quantum_pieri_chains",
    "schubmult_double",
    "schubmult_py",
    "schubmult_q",
    "schubmult_q_double",
    "separated_descents_coeffs",
    "separated_descents_grothmult_double",
    "single_variable_groth",
]
