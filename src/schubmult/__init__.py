__version__ = "3.0.2dev1"
#schubpoly_from_elems?
from .perm_lib import Permutation, permtrim, theta, uncode
from .rings.poly_lib import divide_out_diff, efficient_subs, elem_sym_poly, elem_sym_poly_q, q_vector, split_up, xreplace_genvars
from .rings.quantum_schubert_ring import (
    QDSx,
    QPDSx,
    QPSx,
    QSx,
)
from .rings.schub_poly import div_diff, schubpoly
from .rings.schubert_ring import (
    DSx,
    Sx,
)
from .rings.symmetric_polynomials.complete_sym import CompleteSym
from .rings.symmetric_polynomials.elem_sym import ElemSym
from .rings.symmetric_polynomials.functions import canonicalize_elem_syms, elem_sym_unify, split_out_vars
from .rings.variables import CustomGeneratingSet, GeneratingSet, GeneratingSet_base, MaskedGeneratingSet
from .schub_lib.double import (
    mult_poly_double,
    schub_coprod_double,
    schubmult_double,
    schubmult_double_pair,
    schubmult_double_pair_generic,
)
from .schub_lib.positivity import compute_positive_rep, posify, schubmult_generic_partial_posify
from .schub_lib.quantum import mult_poly_q, schubmult_q, schubmult_q_fast
from .schub_lib.quantum_double import (
    apply_peterson_woodward,
    factor_out_q_keep_factored,
    mult_poly_q_double,
    nil_hecke,
    q_posify,
    schubmult_q_double,
    schubmult_q_double_fast,
    schubmult_q_double_pair_generic,
    schubmult_q_generic_partial_posify,
    schubpoly_quantum,
)
from .schub_lib.schub_lib import check_blocks
from .schub_lib.single import (
    mult_poly_py,
    schub_coprod_py,
    schubmult_py,
)

__all__ = [
    "CompleteSym",
    "CustomGeneratingSet",
    "DSx",
    "ElemSym",
    "GeneratingSet",
    "GeneratingSet_base",
    "MaskedGeneratingSet",
    "Permutation",
    "QDSx",
    "QPDSx",
    "QPSx",
    "QSx",
    "Sx",
    "apply_peterson_woodward",
    "canonicalize_elem_syms",
    "check_blocks",
    "compute_positive_rep",
    "div_diff",
    "divide_out_diff",
    "efficient_subs",
    "elem_sym_poly",
    "elem_sym_poly_q",
    "elem_sym_unify",
    "factor_out_q_keep_factored",
    "mult_poly_double",
    "mult_poly_py",
    "mult_poly_q",
    "mult_poly_q_double",
    "nil_hecke",
    "permtrim",
    "posify",
    "q_posify",
    "q_vector",
    "rings",
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
    "schubmult_q_double_pair_generic",
    "schubmult_q_fast",
    "schubmult_q_generic_partial_posify",
    "schubpoly",
    "schubpoly_quantum",
    "split_out_vars",
    "split_up",
    "theta",
    "uncode",
    "xreplace_genvars",
]
