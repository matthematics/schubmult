from functools import cache

from symengine import sympify

from schubmult.poly_lib.variables import GeneratingSet

NoneVar = 1e10
ZeroVar = 0





@cache
def poly_ring(v: str):
    if v == ZeroVar:
        return tuple([sympify(0) for i in range(100)])
    if v == NoneVar:
        return tuple([sympify(0) for i in range(100)])
    return GeneratingSet(str(v))

# def _schubifyit(func):
#     @wraps(func)
#     def wrapper(f, g):
#         g = _sympify(g)
#         if isinstance(g, Poly):
#             return func(f, g)
#         elif isinstance(g, Integer):
#             g = f.from_expr(g, *f.gens, domain=f.domain)
#             return func(f, g)
#         elif isinstance(g, Expr):
#             try:
#                 g = f.from_expr(g, *f.gens)
#             except PolynomialError:
#                 if g.is_Matrix:
#                     return NotImplemented
#                 expr_method = getattr(f.as_expr(), func.__name__)
#                 result = expr_method(g)
#                 if result is not NotImplemented:
#                     sympy_deprecation_warning(
#                         """
#                         Mixing Poly with non-polynomial expressions in binary
#                         operations is deprecated. Either explicitly convert
#                         the non-Poly operand to a Poly with as_poly() or
#                         convert the Poly to an Expr with as_expr().
#                         """,
#                         deprecated_since_version="1.6",
#                         active_deprecations_target="deprecated-poly-nonpoly-binary-operations",
#                     )
#                 return result
#             else:
#                 return func(f, g)
#         else:
#             return NotImplemented
#     return wrapper
