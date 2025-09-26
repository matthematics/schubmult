# LR rule verification script

import os
import shutil
import sys
import time
from functools import cache
from itertools import zip_longest
from json import dump, load
from math import perm
from multiprocessing import Lock, Manager, Pool, Process, cpu_count

from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.rings.rc_graph_module import FA, RCGraph, try_lr_module, try_lr_module_inject
from schubmult.rings.schubert_ring import Sx

all_rc_graphs = {}

# def generate_rc_graphs(perms, n):
#     global all_rc_graphs
#     for perm in perms:
#         all_rc_graphs[perm] = []
#         for length in range(n):
#             all_rc_graphs[perm].append({g: i for i, g in enumerate(tuple(sorted(RCGraph.all_rc_graphs(perm, length), key=lambda a: a.length_vector())))})


def generate_rc_graph_set(perm, length, bruhat_perm=None):

    from schubmult import FA, ASx, uncode
    from schubmult.rings.rc_graph_module import RCGraph
    if not bruhat_perm:
        bruhat_perm = perm

    if length < len(perm.trimcode):
        return set()
    pants = [*perm.trimcode, *((0,)* (length - len(perm.trimcode)))]
    assert len(pants) == length
    return set((rc[0],rc[1]) for rc in (FA(*pants).coproduct()*(RCGraph()@RCGraph())).value_dict.keys())
    # if len(perm.trimcode) == 0:
    #     return set()
    # lower_elem = ASx(uncode(perm.trimcode[1:]), length - 1)
    # upper_elem = ASx(uncode([perm.trimcode[0]]), 1) * lower_elem
    # return set(key[0] for key in upper_elem.keys() if key[0] != perm)

@cache
def co_principal_perms(perm, length):
    from schubmult import ASx, uncode

    if len(perm.trimcode) == 0:
        return []
    if length < len(perm.trimcode):
        return []
    lower_elem = ASx(uncode(perm.trimcode[1:]), length - 1)
    upper_elem = ASx(uncode([perm.trimcode[0]]), 1) * lower_elem
    return [key[0] for key in upper_elem.keys() if key[0] != perm]

@cache
def code_len(perm, length):
    cd = perm.trimcode
    while len(cd) < length:
        cd = (*cd, 0)
    return tuple(cd)

@cache
def co_principal_seqs(seq):
    from schubmult import ASx, uncode

    if len(seq) == 0:
        return []
    perm = uncode(seq)
    lower_elem = ASx(uncode(seq[1:]), len(seq) - 1)
    upper_elem = ASx(uncode([seq[0]]), 1) * lower_elem
    return [code_len(key[0], len(seq)) for key in upper_elem.keys() if key[0] != perm]

@cache
def good_for_perm(perm, length):
    from schubmult import uncode
    if perm.inv == 0 and length == 0:
        return { (RCGraph(), RCGraph()) }
    if length < len(perm.trimcode):
        raise ValueError(f"Length {length} less than {len(perm.trimcode)} for {perm.trimcode}")
    #print(f"{perm.trimcode,length}")
    
    old_good = good_for_perm(uncode(perm.trimcode[1:]), length - 1)

    new_good = set()
    for (rc1, rc2) in old_good:
        new_good.update([(k[0],k[1]) for k in (FA(perm.trimcode[0]).coproduct() * (rc1@rc2)).value_dict.keys() if k[0].perm.bruhat_leq(perm) and k[1].perm.bruhat_leq(perm)])

    good = set(new_good)
    for (perm2, _) in co_principal_perms(perm, length):
        for (rc1b, rc2b) in good_for_perm(perm2, length):
            good2 = set(good)
            for rc1, rc2 in good2:
                if (rc1.perm == rc1b.perm and rc2.perm == rc2b.perm) and (rc1.length_vector(),rc2.length_vector()) >= (rc1b.length_vector(), rc2b.length_vector()):
                    good = set(good2)
                    good.remove((rc1, rc2))
                    #break

    return good



def build_rcs(seq):
    from schubmult import uncode
    perm = uncode(seq)
    return set((rc[0],rc[1]) for rc in (FA(*seq).coproduct()*(RCGraph()@RCGraph())).value_dict.keys())

def iterative_construct(seq):
    from schubmult import FA
    if len(seq) == 0:
        return {(RCGraph(), RCGraph(), RCGraph())}
    #seqs2 = set(FA(*list(reversed(seq))).coproduct().keys())
    ret_triples = set()
    rev = list(reversed(seq))
    these_triples = {(RCGraph(), RCGraph(), RCGraph())}
    for size in rev:
        these_triples2 = set()
        for rc1, rc2, rc3 in these_triples:
            rc1_new = next(iter(rc1.iterative_act(0, insert=True)))
            rc2_new = next(iter(rc2.iterative_act(0, insert=True)))
            rc3_new = next(iter(rc3.iterative_act(0, insert=True)))
            buildup1 = {(rc1_new, rc2_new, rc3_new)}
            for j in range(size):
                buildup1_new = set()
                for rc11, rc22, rc33 in buildup1:
                    rc1_new_set = rc11.iterative_act(1, insert=False)
                    rc2_new_set = rc22.iterative_act(1, insert=False)
                    rc3_new_set = rc33.iterative_act(1, insert=False)
                    for rc1n in rc1_new_set:
                        for rc3n in rc3_new_set:
                            buildup1_new.add((rc1n, rc22, rc3n))
                    for rc2n in rc2_new_set:
                        for rc3n in rc3_new_set:
                            buildup1_new.add((rc11, rc2n, rc3n))
                buildup1 = buildup1_new
            these_triples2.update(buildup1)
        these_triples = these_triples2
        #print(these_triples)
    ret_triples.update(these_triples)
    return ret_triples



# def rc_less(rc_pair1, rc_pair2):
#     rc1_trns, rc2_trns = rc_pair1[0].transpose(), rc_pair1[1].transpose()
#     rc1b_trns, rc2b_trns = rc_pair2[0].transpose(), rc_pair2[1].transpose()
    # print()
    # print(r)
    # return (rc1_trns.length_vector(), rc2_trns.length_vector()) > (rc1b_trns.length_vector(), rc2b_trns.length_vector())

# @cache
# def binary_good(rc1, rc2):
#     from schubmult import uncode
#     if len(rc1) != len(rc2):
#         return False
#     if len(rc1) == 0 and len(rc2) == 0:
#         return True
#     seq = tuple(a+b for a, b in zip(rc1.length_vector(), rc2.length_vector()))
#     perm = uncode(seq)
#     if not rc1.perm.bruhat_leq(perm) or not rc2.perm.bruhat_leq(perm):
#         return False
#     #print(rc1@rc2)
#     #print("Testing")
#     lower_rc1 = rc1.rowrange(1, len(rc1))
#     lower_rc2 = rc2.rowrange(1, len(rc2))
#     #print(lower_rc1@lower_rc2)
#     if not binary_good(lower_rc1, lower_rc2):
#         #print("Bad")
#         return False
    
#     for seq2 in co_principal_seqs(seq):
#         print(f"For {seq} checking {seq2}")
#         rcs = get_rcs(seq2)
#         for rcb1, rcb2 in rcs:
#             if not binary_good(rcb1, rcb2):
#                 continue
#             if (rc1.perm == rcb1.perm and rc2.perm == rcb2.perm) and rc_less((rc1,rc2),(rcb1,rcb2)):
#             #((rcb1.length_vector(), rcb2.length_vector()) >= (rc1.length_vector(), rc2.length_vector())):
#                 #if
#             #(rc1.length_vector() < rcb1.length_vector() or (rc2.length_vector() < rcb2.length_vector())):
#                 print("Checked")
#                 print(rc1@rc2)
#                 print("Bad becuase")
#                 print(rcb1@rcb2)
#                 print(f"{(rc1.length_vector(),rc2.length_vector())} >= {(rcb1.length_vector(), rcb2.length_vector())}")
#                 return False
#     return True

# def iterative_construct(rc):



def main():
    from schubmult import ASx, Permutation, Sx, uncode
    from schubmult.abc import x, y, z
    from schubmult.perm_lib import Plactic
    from schubmult.rings.rc_graph_module import RCGraph, try_lr_module_biject
    from schubmult.rings.schubert_ring import DoubleSchubertRing
    from schubmult.utils.perm_utils import artin_sequences

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    length = n - 1
    rank = {}
    
    perms.sort(key=lambda p: p.trimcode)
    for perm in perms:
        print("Trying ", perm.trimcode)
        
        #full_coprod = RCGraph.full_rc_coproduct(perm, length)
        full_coprod = RCGraph.principal_rc(perm, length).coproduct()  # Precompute principal rc graphs
        print(full_coprod)

        elem = (ASx@ASx).zero
        for (rc1, rc2), coeff in full_coprod.value_dict.items():
            if uncode([a+b for a,b in zip(rc1.length_vector(), rc2.length_vector())]) == perm:
                elem += coeff * (ASx @ ASx)(((rc1.perm, length), (rc2.perm, length)))
        
        check = ASx(perm, length).coproduct()
        try:
            if perm.inv != 0:
                assert all(v == 0 for v in (elem - check).values())
        except AssertionError:
            print(f"Fail on {perm} at ", time.ctime())
            print(f"{elem=}")
            print(f"{check=}")
            print(f"{(elem - check)=}")
            continue        
        print("Success for ", perm.trimcode)
    # n = int(sys.argv[1])
    # ring = DoubleSchubertRing(z,y)
    # for i0 in range(n):
    #     rc0 = FA(i0).coproduct() * (RCGraph()@RCGraph())
    #     for (rc1, rc2), coeff in rc0.value_dict.items():
    #         print(rc1@rc2)
    #         for perm, coeff in poly_it(rc1.perm, rc2.perm).items():
    #             print(f"  {perm}: {ring(coeff)}")
    #     for i in range(n):
    #         rc = FA(i).coproduct() * (RCGraph()@RCGraph())
            
            

#    rel = generate_co_prinipal_relation(perms, max_length=n)
    
    
    # goodness_seq = {(): {(RCGraph(),RCGraph())}}

    # for length in range(1, n):
    #     this_seq = {}
    #     seqs = artin_sequences(length)
    #     seqs = list(sorted(seqs))
    #     for seq in seqs:
    #         this_seq[seq] = set()
    #         for rc1, rc2 in goodness_seq[seq[1:]]:
    #             this_seq[seq].update({(k[0], k[1]) for k in (FA(seq[0]).coproduct() * (rc1@rc2)).keys() if k[0].perm.bruhat_leq(uncode(seq)) and k[1].perm.bruhat_leq(uncode(seq))})

    #     for seq in seqs:
    #         for (perm2, _) in co_principal_perms(uncode(seq), length):
    #             if len(perm2) > length + 1:
    #                 continue
    #             seq2 = tuple(code_len(perm2, length))
    #             for (rc1b, rc2b) in this_seq[seq2]:
    #                 good2 = set(this_seq[seq])
    #                 for rc1, rc2 in good2:
    #                     if (rc1.perm == rc1b.perm and rc2.perm == rc2b.perm) and (rc1.length_vector() >= rc1b.length_vector() or rc2.length_vector() >= rc2b.length_vector()):
    #                         this_seq[seq] = set(good2)
    #                         this_seq[seq].remove((rc1, rc2))
    #                         break
                    
    #     goodness_seq = this_seq

    # goodness = {uncode(seq): good for seq, good in goodness_seq.items()}

    # for perm in perms:
    #     print(f"Trying {perm.trimcode}")
    #     cop = (ASx@ASx).zero
    #     # good = good_for_perm(perm, length=len(perm.trimcode))
    #     #good = good_for_perm(perm, length=len(perm.trimcode))
    #     rcs = get_rcs(tuple(perm.trimcode))
    #     print(rcs)
    #     sping = iter([rc0[0]@rc0[1] for rc0 in rcs if binary_good(rc0[0], rc0[1])])
    #     dingbat_module = next(sping)
    #     for term in sping:
    #         dingbat_module += term
    #     #dingbat_module = sum([rc0[0]@rc0[1] for rc0 in good.keys()])
    #     print("Correct")
    #     print(try_lr_module(perm, length=len(perm.trimcode)))
    #     print("Test")
    #     print(dingbat_module)
    #     #cop = sum([(ASx@ASx)(((rc0[0].perm,len(rc0[0])),(rc0[1].perm,len(rc0[1])))) for rc0 in good])
    #     cop = dingbat_module.asdtype(ASx@ASx)
    #     test_cop = ASx(perm).coproduct()
    #     tester = test_cop - cop
    
    #             #good = {rc0: g for rc0, g  in at_least_one.items() if g}
    #     # print("Correct")
    #     # print(try_lr_module(perm, length=len(perm.trimcode)))
    #     # print("Test")
        
    #     try:
    #         assert all(v == 0 for v in tester.values()), f"Failed for {perm=} {cop=} {test_cop=}"
    #     except AssertionError as e:
    #         print(f"{cop=}")
    #         print(f"{test_cop=}")
    #         print("Difference:")
    #         for k, v in tester.items():
    #             if v != 0:
    #                 print(f"  {k}: {v}")
    #         print("Failure for ", perm.trimcode)
    #         raise
    #     print("Success for ", perm.trimcode)
        #goodness[perm] = good
    # try_sum = {}
    # for perm in perms:
    #     try_sum[perm] = try_sum.get(perm, {})
    #     for (rc1, rc2) in goodness[perm]:
    #         try_sum[perm][(rc1.perm, rc2.perm)] = try_sum.get((rc1.perm, rc2.perm), 0) + rc1.polyvalue(x) * rc2.polyvalue(x)
    #     rc_pairs = generate_rc_graph_set(perm, length=len(perm.trimcode))
    #     for (perm2, _) in co_principal_perms(perm, len(perm.trimcode)):
    #         length = len(perm.trimcode)
    #         for (rc1b, rc2b) in good_for_perm(perm2, length):
    #             for rc1, rc2 in rc_pairs:
    #                 if (rc1.perm == rc1b.perm and rc2.perm == rc2b.perm) and (rc1.length_vector() >= rc1b.length_vector() or rc2.length_vector() >= rc2b.length_vector()):
    #                     try_sum[perm2]=try_sum.get(perm2,{})
    #                     try_sum[perm2][(rc1.perm, rc2.perm)] = try_sum[perm2].get((rc1.perm, rc2.perm), 0) + rc1.polyvalue(x) * rc2.polyvalue(x)
    # Permutation.print_as_code=True
    # for perm in perms:
    #     for (perm1, perm2), coeff in try_sum[perm].items():
    #         print(f"{perm.trimcode,perm1.trimcode,perm2.trimcode}: {Sx(coeff)}")
    #         assert Sx(coeff).get(perm, 0) == (Sx(perm1) * Sx(perm2)).get(perm, 0), f"Failed {perm.trimcode,perm1.trimcode,perm2.trimcode}: {Sx(coeff)} vs {Sx(perm1) * Sx(perm2)}"
if __name__ == "__main__":
    main()
