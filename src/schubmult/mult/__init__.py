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
Exports resolve lazily (PEP 562) so that importing one kernel does not load all the others.
"""

_exports = {
    "mult_poly_double": ("double", "mult_poly_double"),
    "schubmult_double": ("double", "schubmult_double"),
    "grothmult_py": ("groth", "grothmult_py"),
    "mult_poly_groth": ("groth", "mult_poly_groth"),
    "grothmult_q": ("groth_quantum", "grothmult_q"),
    "grothmult_q_dict": ("groth_quantum", "grothmult_q_dict"),
    "grothmult_q_pieri": ("groth_quantum", "grothmult_q_pieri"),
    "compute_positive_rep": ("positivity", "compute_positive_rep"),
    "posify": ("positivity", "posify"),
    "schubmult_q": ("quantum", "schubmult_q"),
    "factor_out_q": ("quantum_double", "factor_out_q"),
    "schubmult_q_double": ("quantum_double", "schubmult_q_double"),
    "separated_descents_grothmult_double": ("separated_descents", "grothmult_double"),
    "separated_descents_coeffs": ("separated_descents", "separated_descents_coeffs"),
    "mult_poly_py": ("single", "mult_poly_py"),
    "schubmult_py": ("single", "schubmult_py"),
}
_exports.update(
    {
        name: ("groth_double", name)
        for name in (
            "elem_sym_perms_groth",
            "epsilon_chain",
            "groth_elem_sym_poly",
            "grothmult_double",
            "grothmult_double_block",
            "grothmult_double_pieri",
            "monk_chain",
            "mult_poly_groth_double",
            "one_plus_beta_x_groth",
            "single_variable_groth",
        )
    },
)
_exports.update(
    {
        name: ("groth_quantum_double", name)
        for name in (
            "groth_elem_sym_poly_q",
            "grothmult_q_double",
            "grothmult_q_double_dict",
            "grothmult_q_double_pieri",
            "grothmult_q_double_top",
            "lm_quantize",
            "qgroth_poly",
            "quantum_elem_sym",
            "quantum_pieri_chains",
        )
    },
)


def __getattr__(name):
    if name not in _exports:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    modname, attr = _exports[name]
    val = getattr(importlib.import_module(f"{__name__}.{modname}"), attr)
    globals()[name] = val
    return val


def __dir__():
    return sorted(set(globals()) | set(__all__))


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
