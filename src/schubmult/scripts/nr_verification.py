# LR rule verification script

import os
import shutil
import sys
import time
from functools import cache
from itertools import zip_longest
from json import dump, load
from multiprocessing import Lock, Manager, Pool, Process, cpu_count

from schubmult.rings.free_algebra_basis import WordBasis
from schubmult.rings.rc_graph_module import FA, RCGraph, try_lr_module
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
    return [key for key in upper_elem.keys() if key[0] != perm]

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
                if (rc1.perm == rc1b.perm and rc2.perm == rc2b.perm) and (rc1.length_vector() >= rc1b.length_vector() or rc2.length_vector() >= rc2b.length_vector()):
                    good = set(good2)
                    good.remove((rc1, rc2))
                    #break

    return good

def code_len(perm, length):
    cd = perm.trimcode
    while len(cd) < length:
        cd = (*cd, 0)
    return cd

def main():
    from schubmult import ASx, Permutation, Sx, uncode
    from schubmult.abc import x
    from schubmult.rings.rc_graph_module import RCGraph
    from schubmult.utils.perm_utils import artin_sequences

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)

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

    for perm in perms:
        print(f"Trying {perm.trimcode}")
        cop = (ASx@ASx).zero
        # good = good_for_perm(perm, length=len(perm.trimcode))
        good = good_for_perm(perm, length=len(perm.trimcode))
        sping = iter([rc0[0]@rc0[1] for rc0 in good])
        dingbat_module = next(sping)
        for term in sping:
            dingbat_module += term
        #dingbat_module = sum([rc0[0]@rc0[1] for rc0 in good.keys()])
        print("Correct")
        print(try_lr_module(perm, length=len(perm.trimcode)))
        print("Test")
        print(dingbat_module)
        #cop = sum([(ASx@ASx)(((rc0[0].perm,len(rc0[0])),(rc0[1].perm,len(rc0[1])))) for rc0 in good])
        cop = dingbat_module.asdtype(ASx@ASx)
        test_cop = ASx(perm).coproduct()
        tester = test_cop - cop
    
                #good = {rc0: g for rc0, g  in at_least_one.items() if g}
        # print("Correct")
        # print(try_lr_module(perm, length=len(perm.trimcode)))
        # print("Test")
        
        try:
            assert all(v == 0 for v in tester.values()), f"Failed for {perm=} {cop=} {test_cop=}"
        except AssertionError as e:
            print(f"{cop=}")
            print(f"{test_cop=}")
            print("Difference:")
            for k, v in tester.items():
                if v != 0:
                    print(f"  {k}: {v}")
            print("Failure for ", perm.trimcode)
            raise
        print("Success for ", perm.trimcode)
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
