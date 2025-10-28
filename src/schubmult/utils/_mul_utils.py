from schubmult.symbolic import sympify
from schubmult.utils.perm_utils import add_perm_dict


def _mul_schub_dicts(dict1, dict2, basis1, basis2, best_effort_positive=False):
    this_dict = {}
    for k, v in dict2.items():
        for kd, vd in dict1.items():
            did_positive = False
            to_mul = v * vd
            if best_effort_positive:
                try:
                    this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_positive_product(kd, k, basis2).items()})
                    did_positive = True
                except Exception:
                    raise
            if not did_positive:
                this_dict = add_perm_dict(this_dict, {k1: v1 * to_mul for k1, v1 in basis1.cached_product(kd, k, basis2).items()})
    return this_dict


def _tensor_product_of_dicts_first(d1, d2):
    ret_dict = {}
    for k1, v1 in d1.items():
        this_dict = {}
        v1 = sympify(v1)
        for k2, v2 in d2.items():
            v2 = sympify(v2)
            this_dict[(k1, k2)] = v1 * v2
        ret_dict = add_perm_dict(ret_dict, this_dict)
    return ret_dict


def _tensor_product_of_dicts(d1, d2):
    ret_dict = {}
    for k1, v1 in d1.items():
        this_dict = {}
        v1 = sympify(v1)
        for k2, v2 in d2.items():
            v2 = sympify(v2)
            if isinstance(k1, tuple):
                this_dict[(*k1, k2)] = v1 * v2
            else:
                this_dict[(k1, k2)] = v1 * v2
        ret_dict = add_perm_dict(ret_dict, this_dict)
    return ret_dict
