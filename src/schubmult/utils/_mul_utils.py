"""Dict-level helpers for ring multiplication and tensor products of ``{key: coeff}`` expansions."""

from schubmult.symbolic import sympify
from schubmult.utils.perm_utils import add_perm_dict


def _mul_schub_dicts(dict1, dict2, basis1, basis2, best_effort_positive=False):
    """Bilinear product of two expansions using ``basis1.cached_product(k1, k2, basis2)`` (or the
    ``cached_positive_product`` variant) for each pair of keys.
    """
    this_dict = {}
    for k, v in dict2.items():
        for kd, vd in dict1.items():
            to_mul = v * vd
            prod = basis1.cached_positive_product(kd, k, basis2) if best_effort_positive else basis1.cached_product(kd, k, basis2)
            if to_mul == 1:
                for k1, v1 in prod.items():
                    this_dict[k1] = this_dict[k1] + v1 if k1 in this_dict else v1
            else:
                for k1, v1 in prod.items():
                    this_dict[k1] = this_dict[k1] + v1 * to_mul if k1 in this_dict else v1 * to_mul
    return this_dict


def _tensor_product_of_dicts_first(d1, d2):
    """``{(k1, k2): v1 * v2}`` over all pairs; keys are always wrapped as a pair."""
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
    """Like `_tensor_product_of_dicts_first` but flattens: a tuple key ``k1`` becomes ``(*k1, k2)``."""
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
