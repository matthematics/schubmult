from schubmult import *
from schubmult.utils.schub_lib import compute_vpathdicts
from schubmult.utils.perm_utils import add_perm_dict


def count_paths(perm, length):
    numvars = length
    res = {}
    def word_elem(p, k, *args):
        return FA(k - p)
    expr = {k[1:]: v for k, v in Sx(perm * ~uncode(list(range(perm.inv + numvars, perm.inv, -1)))).in_SEM_basis(elem_func=word_elem).items()}
    return expr

if __name__ == "__main__":
    import sys
    from schubmult.utils.schub_lib import complete_sym_positional_perms_down
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        if perm.inv == 0:
            continue
        for length in range(len(perm.trimcode), n):
            # mulperm = perm.mul_dominant()
            # code_to_work = mulperm.trimcode
            # path_length = len(code_to_work)
            # assert path_length <= length, f"Path length {path_length} is greater than length {length} for permutation {perm}"
            # if path_length < length:
            #     while path_length < length:
            #         code_to_work += [code_to_work[-1]]
            #         path_length += 1
            # stack = [[[(1, mulperm)], [(), 1]]]
            # the_coeffs = {}
            # while stack:
            #     the_path, (weight, path_sign) = stack.pop()
                
            #     if the_path[0][0] == 0 and the_path[0][1] == ~perm:
            #         the_coeffs[tuple(weight)] = the_coeffs.get(tuple(weight), 0) + path_sign
            #     elif the_path[0][0] > 0:
            #         for next_perm, weight_diff, the_sign in kdown_perms(the_path[0][1], code_to_work[the_path[0][0]], the_path[0][0]):
            #             if (~perm).bruhat_leq(next_perm):
            #                 stack.append([[(the_path[0][0] + 1, next_perm)] + the_path, [(weight_diff,) + weight, path_sign * the_sign]])
            the_coeffs = count_paths(perm, length)
            dual_schub = ASx(perm, length).change_basis(WordBasis)
            for k, v in the_coeffs.items():
                try:
                    assert dual_schub.get(k, 0) == v, f"Failed for {perm} at length {length} with weight {k} and coeff {v}, got {dual_schub.get(k, 0)}, {the_coeffs=} {dual_schub=}"
                except AssertionError as e:
                    print(e)
                    print(f"Permutation: {perm}, length: {length}, weight: {k}, coeff: {v}, dual_schub coeff: {dual_schub.get(k, 0)}")
                    #raise