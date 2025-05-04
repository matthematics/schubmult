from functools import cache

import sympy
from symengine import SympifyError, sympify

from schubmult.rings.variables import CustomGeneratingSet, GeneratingSet
from schubmult.utils.perm_utils import add_perm_dict

NoneVar = 1e10
ZeroVar = 0


class NotEnoughGeneratorsError(ValueError):
    pass


@cache
def poly_ring(v: str):
    if v == ZeroVar:
        return CustomGeneratingSet(tuple([sympify(0) for i in range(100)]))
    if v == NoneVar:
        return CustomGeneratingSet(tuple([sympify(0) for i in range(100)]))
    return GeneratingSet(str(v))


def _mul_schub_dicts(dict1, dict2, basis1, basis2, best_effort_positive=True, _sympify=False):
    this_dict = {}
    symp = lambda x: x  # noqa: E731
    # if _sympify:
    #     symp = sympy.sympify
    for k, v in dict2.items():
        for kd, vd in dict1.items():
            did_positive = False
            to_mul = symp(v) * symp(vd)
            if best_effort_positive:
                try:
                    this_dict = add_perm_dict(this_dict, {k1: symp(v1) * to_mul for k1, v1 in basis1.cached_positive_product(kd, k, basis2).items()})
                    did_positive = True
                except Exception:
                    did_positive = False
            if not did_positive:
                this_dict = add_perm_dict(this_dict, {k1: symp(v1) * to_mul for k1, v1 in basis1.cached_product(kd, k, basis2).items()})
    return this_dict


def _tensor_product_of_dicts(d1, d2):
    ret_dict = {}
    for k1, v1 in d1.items():
        this_dict = {}
        try:
            v1 = sympify(v1)
        except SympifyError:
            v1 = sympy.sympify(v1)
        for k2, v2 in d2.items():
                try:
                    v2 = sympify(v2)
                except SympifyError:
                    v2 = sympy.sympify(v2)
                    v1 = sympy.sympify(v1)
                if isinstance(k1, tuple):
                    this_dict[(*k1, k2)] = v1 * v2
                else:
                    this_dict[(k1, k2)] = v1 * v2
        ret_dict = add_perm_dict(ret_dict, this_dict)
    return ret_dict
