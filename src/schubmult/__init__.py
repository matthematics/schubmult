__version__ = "2.1.0"

from .poly_lib import GeneratingSet, base_index, div_diff, efficient_subs, elem_sym_poly, elem_sym_poly_q, q_vector, schubpoly, skew_div_diff, xreplace_genvars
from .schub_lib.double import (
    compute_positive_rep,
    mult_poly_double,
    schub_coprod_double,
    schubmult_double,
    schubmult_double_pair,
    schubmult_double_pair_generic,
    schubmult_generic_partial_posify,
)
from .schub_lib.quantum import mult_poly_q, schubmult_q, schubmult_q_fast
from .schub_lib.quantum_double import factor_out_q_keep_factored, mult_poly_q_double, schubmult_q_double, schubmult_q_double_fast, schubpoly_quantum
from .schub_lib.single import (
    mult_poly_py,
    schub_coprod_py,
    schubmult_py,
)

__all__ = [
    "GeneratingSet",
    "base_index",
    "compute_positive_rep",
    "mult_poly_double",
    "schub_coprod_double",
    "schubmult_double",
    "schubmult_double_pair",
    "schubmult_double_pair_generic",
    "schubmult_generic_partial_posify",
    "mult_poly_py",
    "schub_coprod_py",
    "schubmult_py",
    "mult_poly_q",
    "schubmult_q",
    "schubmult_q_fast",
    "factor_out_q_keep_factored",
    "mult_poly_q_double",
    "schubmult_q_double",
    "schubmult_q_double_fast",
    "schubpoly_quantum",
    "div_diff",
    "efficient_subs",
    "elem_sym_poly",
    "elem_sym_poly_q",
    "perm_act",
    "q_vector",
    "schubpoly",
    "skew_div_diff",
    "xreplace_genvars",
]
