from sage.combinat.composition import Composition
from sage.combinat.permutation import Permutation

from schubmult.perm_lib import permtrim, trimcode, uncode


def _coerce_index(indexed_obj, is_comp, should_be_comp):
    if is_comp == should_be_comp:
        if isinstance(indexed_obj, list) or isinstance(indexed_obj, tuple):
            if is_comp:
                return Composition(trimcode(permtrim(uncode(list(indexed_obj)))))
            return Permutation(permtrim(list(indexed_obj))).remove_extra_fixed_points()
        if not is_comp:
            if isinstance(indexed_obj, Permutation):
                return indexed_obj.remove_extra_fixed_points()
            if isinstance(indexed_obj, dict):
                {
                Permutation(permtrim(list(k))).remove_extra_fixed_points(): v
                for k, v in indexed_obj.items()
            }
    elif is_comp:
        if (
            isinstance(indexed_obj, list)
            or isinstance(indexed_obj, tuple)
            or isinstance(indexed_obj, Composition)
        ):
            return Permutation(permtrim(uncode(list(indexed_obj)))).remove_extra_fixed_points()

        if isinstance(indexed_obj, dict):  # keys are comps
            return {
                Permutation(permtrim(uncode(list(k)))).remove_extra_fixed_points(): v
                for k, v in indexed_obj.items()
            }
    else:
        if (
            isinstance(indexed_obj, list)
            or isinstance(indexed_obj, tuple)
            or isinstance(indexed_obj, Permutation)
        ):
            return Composition(trimcode(list(indexed_obj)))

        if isinstance(indexed_obj, dict):  # keys are comps
            return {
                Composition(trimcode(permtrim(list(k)))): v
                for k, v in indexed_obj.items()
            }
    raise TypeError
