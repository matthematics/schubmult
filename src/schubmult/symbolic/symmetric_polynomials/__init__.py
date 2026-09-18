
"""Symbolic (factorial) elementary and complete symmetric polynomials as SymPy functions.

``E(p, k, *vars)``/``e`` and ``H(p, k, *vars)``/``h`` are the elementary and complete symmetric
polynomials of degree ``p`` in the first ``k`` generators, kept unevaluated so Schubert
expansions can be written in the SEM basis; the ``Factorial*`` variants carry a second
variable set. `functions` holds canonicalization and variable-splitting utilities, and
`qelem_sym` the quantum elementary symmetric polynomials.
"""

from .complete_sym import CompleteSym, CompleteSym_base, FactorialCompleteSym, H, h
from .elem_sym import E, ElemSym, ElemSym_base, FactorialElemSym, e
from .functions import canonicalize_elem_syms, canonicalize_elem_syms_coeff, coeffvars, degree, genvars, is_of_func_type, numvars, split_out_vars
from .qelem_sym import E_q, QFactorialElemSym
