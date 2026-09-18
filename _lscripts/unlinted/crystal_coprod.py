from time import time

from matplotlib.pyplot import rc
from schubmult import RCGraph
from schubmult import RCGraphRing, ring_elem_to_highest_weight, tensor_to_highest_weight, tensor_to_highest_weight2
from sympy import pretty_print

from schubmult import *


def try_lr_module0(perm, length=None):
    from schubmult import RCGraphRing

    from schubmult import ASx, uncode
    ring = RCGraphRing()
    tring = ring @ ring
    # print(f"Starting {perm}")
    if length is None:
        length = len(perm.trimcode)
    elif length < len(perm.trimcode):
        raise ValueError("Length too short")
    if perm.inv == 0:
        return tring((RCGraph([()]*length),RCGraph([()]*length)))
    lower_perm = uncode(perm.trimcode[1:])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = try_lr_module(lower_perm, length - 1)
    # assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
    #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
    #  #  # print("Going for it")
    #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
    #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")

    cprod = tring.zero

    for j in range(perm.trimcode[0] + 1):
        cprod += tring.ext_multiply(ring(RCGraph.one_row(j)), ring(RCGraph.one_row(perm.trimcode[0] - j)))
    #print(cprod)
    ret_elem = cprod * lower_module1
    #  #  # print(f"{ret_elem=}")
    # assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = tensor_to_highest_weight(tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)}))

    if length == 1:
        return ret_elem
    #keys = set(ret_elem.keys())
    # print(f"{repr(keys)=} {perm=}")
    up_elem = ASx(uncode([perm.trimcode[0]]),1) * elem
    # print(f"{up_elem=}")
    for key, coeff in up_elem.items():
        if key[0] != perm:
            assert coeff == 1
            ret_elem -= tensor_to_highest_weight(try_lr_module(key[0], length))
            #    ret_elem -= cff2 * tring(RCGraph.to_highest_weight_pair(rc1_bad, rc2_bad)[0])
                #        break
    # print(f"Done {perm}")
    #ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k in keys})
    # assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})
    return ret_elem

def try_lr_module(perm, length=None):
    from schubmult import RCGraphRing

    from schubmult import ASx, uncode
    ring = RCGraphRing()
    tring = ring @ ring
    # print(f"Starting {perm}")
    if length is None:
        length = len(perm.trimcode)
    elif length < len(perm.trimcode):
        raise ValueError("Length too short")
    if perm.inv == 0:
        return tring((RCGraph([()]*length),RCGraph([()]*length)))
    if length > len(perm.trimcode):
        mul_elem = 0
        lower_perm = perm
    else:
        mul_elem = perm.trimcode[-1]
        lower_perm = uncode(perm.trimcode[:-1])
    elem = ASx(lower_perm, length - 1)
    lower_module1 = try_lr_module(lower_perm, length - 1)
    # assert isinstance(lower_module1, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"
    #  #  # print(f"Coproducting {ASx(uncode([perm.trimcode[0]]), 1).coproduct()=}")
    #  #  # print(ASx(uncode([perm.trimcode[0]]), 1).coproduct())
    #  #  # print("Going for it")
    #  #  # print(f"{type(lower_module1)=} {lower_module1=}")
    #  #  # print(f"{type(ASx(uncode([perm.trimcode[0]]), 1).coproduct())=}")

    cprod = tring.zero

    for j in range(mul_elem + 1):
        cprod += tring.ext_multiply(ring(RCGraph.one_row(j)), ring(RCGraph.one_row(mul_elem - j)))
    #print(cprod)
    ret_elem = lower_module1 *cprod
    #  #  # print(f"{ret_elem=}")
    # assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(lower_module1)} {lower_perm=} {length=}"

    ret_elem = tensor_to_highest_weight2(tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)}))

    if length == 1:
        return ret_elem
    #keys = set(ret_elem.keys())
    # print(f"{repr(keys)=} {perm=}")
    up_elem = elem * ASx(uncode([mul_elem]),1)
    # print(f"{up_elem=}")
    for key, coeff in up_elem.items():
        if key[0] != perm:
            assert coeff == 1
            ret_elem -= tensor_to_highest_weight2(try_lr_module(key[0], length))
            #    ret_elem -= cff2 * tring(RCGraph.to_highest_weight_pair(rc1_bad, rc2_bad)[0])
                #        break
    # print(f"Done {perm}")
    #ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k in keys})
    # assert isinstance(ret_elem, TensorModule), f"Not TensorModule {type(ret_elem)} {perm.trimcode=}"
    ret_elem = tring.from_dict({k: v for k, v in ret_elem.items() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)})
    return ret_elem


if __name__ == "__main__":
    import sys
    ring = RCGraphRing()
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    tring = ring @ ring
    for perm in perms:
        for length in range(len(perm.trimcode), n):
        # stink_elem = ring_elem_to_highest_weight(ring.from_free_algebra_element(ASx(perm).change_basis(WordBasis)))
        # cprd = ASx(perm).change_basis(WordBasis).coproduct()
        # elem = tring.zero
        # for (w1, w2), coeff in cprd.items():
        #     elem += coeff*tring.ext_multiply(ring.from_free_algebra_element(FA(*w1)),ring.from_free_algebra_element(FA(*w2)))
        # elem = tensor_to_highest_weight(elem)
        # elem = elem.ring.from_dict({(k1, k2): v for (k1, k2), v in elem.items() if k1.perm.bruhat_leq(perm) and k2.perm.bruhat_leq(perm)})
            print(f"Coprod of {perm.trimcode=}")
            try_mod = try_lr_module(perm, length)
            pretty_print(try_mod)
            elem = 0
            for (rc1, rc2), conch_shell in try_mod.items():
                # elem += (rc1 @ rc2).asdtype(ASx @ ASx)
                # print(f"FYI {perm.trimcode} 1")
                # print(rc1)
                # print(f"FYI {perm.trimcode} 2")
                # print(rc2)
                elem += conch_shell*(ASx @ ASx)(((rc1.perm, len(rc1)), (rc2.perm, len(rc2))))
            check = ASx(perm, length).coproduct()
            try:
                assert all(v == 0 for v in (elem - check).values())
            except AssertionError:
                print(f"Fail on {perm}")
                print(f"{elem=}")
                print(f"{check=}")
                print(f"{(elem - check)=}")
                continue
            # print("Stink_elem:")
            # pretty_print(stink_elem)
            print(f"Stinkcess {perm, length}")
