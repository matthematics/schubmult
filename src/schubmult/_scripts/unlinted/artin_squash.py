from schubmult import *
from sympy import pretty_print, S
from schubmult.utils.perm_utils import artin_sequences

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    
    # for seq1, seq2 in itertools.product(artin_sequences(n), repeat=2):
    #     rc1 = r.one
    #     rc2 = r.one
    #     full_rc = r.one
    #     for k in range(1, n):
    #         rc1 = rc1.resize(k)
    #         rc2 = rc2.resize(k)
    #         full_rc = full_rc.resize(k)
    #         new_rc1 = r.zero
    #         new_rc2 = r.zero
    #         new_full_rc = r.zero
    #         for elem_rc1 in RCGraph.all_rc_graphs(uncode([0] * (k - seq1[k - 1]) + [1] * seq1[k - 1]), k):
    #             new_rc1 += r(rc1) % r(elem_rc1)
    #             new_full_rc += r(full_rc) % r(elem_rc1)
    #             for elem_rc2 in RCGraph.all_rc_graphs(uncode([0] * (k - seq2[k - 1]) + [1] * seq2[k - 1]), k):
    #                 new_rc2 += r(rc2) % r(elem_rc2)
    #                 new_full_rc += r(full_rc) % r(elem_rc2)
    #         rc1 = new_rc1
    #         rc2 = new_rc2
    #         full_rc = new_full_rc
    complete_e = [RCGraph.all_rc_graphs(uncode([0] * (n - k) + [1] * k), n) for k in range(n + 1)]

        # rc1 = RCGraph.from_artin_sequence(seq1)
        # rc2 = RCGraph.from_artin_sequence(seq2)
        # new_rc = rc1.squash_product(rc2)
        # print("Seq1:", seq1)
        # print("Seq2:", seq2)
        # print("New RC:")
        # pretty_print(new_rc)