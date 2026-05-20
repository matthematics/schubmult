#!/usr/bin/env python3
"""
Script: make_principal.py
Description: Entry point for principal object construction in schubmult.
"""
import sys
from schubmult import *

def main():
    if len(sys.argv) < 2:
        print("Usage: make_principal <n>")
        sys.exit(1)
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = BoundedRCFactorAlgebra()
    for perm in perms:
        if perm.inv == 0:
            continue
        cd = [*perm.trimcode] + [0] * (n - len(perm.trimcode))
        inv = sum(cd)
        index = n
        spum = []
        while inv > 0:
            #if len(uncode(cd)) < index:
            while index > 0 and cd[index - 1] == 0:
                index -= 1
                #continue
            if index == 0:
                raise ValueError(f"Ran out of index for permutation {perm} with code {cd} and inversion number {inv}")
            #if index != n - 1:
            index -= 1
            spots = [0] * (index)
            
            for i in range(index):
                if cd[i] > 0:
                    spots[i] = 1
                    cd[i] -= 1
            num_spots = sum(spots)        
            inv -= num_spots

            elem_rc = next(iter(RCGraph.all_rc_graphs(uncode([0] * (index - num_spots) + [1] * num_spots), index, weight=tuple(spots))))
            #buildup_rc = buildup_rc.resize(len(elem_rc)).squash_product(elem_rc)
            spum = [elem_rc, *spum]
            
        rc = r.key_to_rc_graph(r.make_key(tuple(spum), n))
        print(f"{repr(rc)}")
        assert rc.perm == perm, f"Principal RC graph does not match original permutation: {rc.perm} != {perm}"
        assert rc.is_principal, f"RC graph is not principal: {rc}"
        #return tuple(spum)
        # tensor = CrystalGraphTensor(*RCGraph.principal_rc_factorization(perm))
        # for rc in [rcc for rcc in RCGraph.all_rc_graphs(perm) if rcc.snap_qy().length_vector == cd]:
        #     rc_high, seq = rc.to_highest_weight()
            
        #     tensor2 = CrystalGraphTensor(*tensor).to_highest_weight()[0].reverse_raise_seq(seq)
            
        #     assert r.key_to_rc_graph(r.make_key(tensor2, len(rc))).resize(len(rc)) == rc, f"Principal RC graph does not match original RC graph {r.key_to_rc_graph(r.make_key(tensor2, len(rc_bot))).resize(len(rc))} != {rc}"
        print("Processed permutation:", perm)

if __name__ == "__main__":
    main()
