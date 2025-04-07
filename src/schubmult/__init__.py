__version__ = "3.0.dev"

from .perm_lib.perm_lib import Permutation
from .poly_lib.poly_lib import efficient_subs, elem_sym_poly, elem_sym_poly_q, q_vector, xreplace_genvars
from .poly_lib.schub_poly import div_diff, schubpoly, skew_div_diff
from .poly_lib.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from .rings._quantum_schubert_polynomial_ring import (
    QDSx,
    QSx,
)
from .rings._schubert_polynomial_ring import (
    DSx,
    Sx,
)
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
    "DSx",
    "GeneratingSet",
    "MaskedGeneratingSet",
    "Sx",
    "compute_positive_rep",
    "CustomGeneratingSet",
    "div_diff",
    "efficient_subs",
    "elem_sym_poly",
    "elem_sym_poly_q",
    "factor_out_q_keep_factored",
    "GeneratingSet_base", 
    "mult_poly_double",
    "mult_poly_py",
    "mult_poly_q",
    "mult_poly_q_double",
    "Permutation",
    "q_vector",
    "QDSx",
    "QSx",
    "schub_coprod_double",
    "schub_coprod_py",
    "schubmult_double",
    "schubmult_double_pair",
    "schubmult_double_pair_generic",
    "schubmult_generic_partial_posify",
    "schubmult_py",
    "schubmult_q",
    "schubmult_q_double",
    "schubmult_q_double_fast",
    "schubmult_q_fast",
    "schubpoly",
    "schubpoly_quantum",
    "skew_div_diff",
    "xreplace_genvars",
]
