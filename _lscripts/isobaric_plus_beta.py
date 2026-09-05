from schubmult import *
from schubmult.symbolic import *
from schubmult.rings.schubert.nil_hecke import NilHeckeRing
from schubmult.abc import x, y
from schubmult.utils.schub_lib import elem_sym_perms_op

nh = NilHeckeRing(y)
nhx = NilHeckeRing(x)

subs_dict = {x[i]: y[i] for i in range(50)}

def _weight_of_perm(perm):
    if perm.inv == 0:
        return S.One
    inv_code = (~perm).trimcode
    result = S.One
    for index, value in enumerate(inv_code, start=1):
        if value != 0:
            # for i in range(value):
            #     result *= (S.One + (Gx._beta) * y[index + i])
            result *= (S.One + (Gx._beta) * y[index]) ** value
    return result

def positive_isobaric(perm):
    #word = perm.inverse_code_word
    word = tuple(reversed((~perm).code_word))
    result = nh.zero
    subwords = Permutation.all_subwords(word)
    for new_word in subwords:
        new_perm = Permutation.hecke_ref_product(*new_word)
        result += ((-Gx._beta)**(len(word) - new_perm.inv) * _weight_of_perm(new_perm))*nh(new_perm) 
    return result

def real_isobaric(perm):
    if perm.inv == 0:
        return nh.one
    d = max(perm.descents()) + 1
    interim_result = nhx.isobaric(perm, groth=True, groth_beta=Gx._beta)
    result = nh.zero
    for perm2, coeff in interim_result.items():
        result += coeff.subs(subs_dict) * nh(perm2)
    return result
    # return interim_result
# def commute_past(perm, genset, elem_start, elem_length, beta):
#     if perm.inv == 0:
#         return prod([(1+beta*genset[i]) for i in range(elem_start, elem_start + elem_length)])*nh.one
#         #return tuple(range(elem_start, elem_start + elem_length)), nh.one
#     d = min(perm.descents()) + 1
#     if d < elem_start:
#         raise ValueError(f"commute_past: perm={perm} elem_start={elem_start} elem_length={elem_length} beta={beta}")
#     result = nh.zero
#     perm_list = elem_sym_perms_op(perm, elem_length, elem_length)
#     for perm_down, _ in perm_list:
#         result += prod([(1+beta * genset[perm[i]]) for i in range(elem_start, elem_start + elem_length) if perm[i] == perm_down[i]]) * nh(perm_down)
#     return result

# def decomp_into_start_length(perm):
#     if perm.inv == 0:
#         return ()
#     index = min([i for i in range((~perm).max_descent) if (~perm).trimcode[i] != 0])
#     d = (~perm).trimcode[index]
#     downperm_inv = uncode((~perm).trimcode[:index] + [0] + (~perm).trimcode[index + 1:])
#     if downperm_inv.inv != 0:
#         return  decomp_into_start_length(downperm_inv) + ((index + 1, d),)
#     return ((index + 1, d),)

# def isobaric_plus_beta(perm, genset, beta):
#     if perm.inv == 0:
#         return nh.one
#     decomp = decomp_into_start_length(perm)
#     print(perm.trimcode)
#     print(decomp)
#     result = nh.one
#     code_buildup = [0] * (len(perm) - 1)
#     for start, length in decomp:
#         result2 = nh.zero
#         for perm2, coeff in result.items():
#             result2 += coeff * commute_past(perm2, genset, start, length, beta)
#         result = nh.zero
#         make_perm = ~uncode(code_buildup)
#         for perm3, coeff in result2.items():
#             if (make_perm * perm3).inv == make_perm.inv+perm3.inv:
#                 result += coeff * nh(make_perm * perm3)
#         print(f"{code_buildup=} {start=}")
#         code_buildup[start - 1] = length
#         result = result2
#     return result

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    beta = Gx._beta
    perms = Permutation.all_permutations(n)
    for perm in perms:
        if perm.inv == 0:
            continue
        pos_test = positive_isobaric(perm)
        real_test = real_isobaric(perm)
        test_test = pos_test - real_test
        test_test = nh.from_dict({k: v.expand() for k, v in test_test.items() if v.expand() != 0})
        print(f"{perm.trimcode}: {test_test}")