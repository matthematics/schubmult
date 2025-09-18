# LR rule verification script

import os
import shutil
import sys
import time
from json import dump, load
from multiprocessing import Lock, Manager, Pool, Process, cpu_count

from schubmult.rings.rc_graph_module import RCGraph, try_lr_module


def co_principal_perms(perm):
    from schubmult import ASx, uncode

    if len(perm.trimcode) == 0:
        return []
    lower_elem = ASx(uncode(perm.trimcode[1:]))
    upper_elem = ASx(uncode([perm.trimcode[0]]), 1) * lower_elem
    return [key[0] for key in upper_elem.keys() if key[0] != perm]

def main():
    from schubmult import Permutation

    # try:
    all_rc_graphs = {}

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    perms.sort(key=lambda p: p.trimcode, reverse=True)
    latest_seen_left = {perm: [{} for _ in range(n)] for perm in perms}
    latest_seen_right = {perm: [{} for _ in range(n)] for perm in perms}
    for perm in perms:
        all_rc_graphs[perm] = []
        for length in range(n):
            all_rc_graphs[perm].append({g: i for i, g in enumerate(tuple(sorted(RCGraph.all_rc_graphs(perm, length), key=lambda a: a.length_vector())))})
    for perm in perms:
        modules = try_lr_module(perm)
        print(all_rc_graphs[perm])
        for (rc1, rc2), coeff in modules.items():
            print(f"{perm.trimcode}: {rc1.perm.trimcode} {len(rc2.perm.trimcode)}")
            print(f"  {latest_seen_left[perm][len(rc1)].get(rc1.perm, -1)} * {latest_seen_right[perm][len(rc2)].get(rc2.perm, -1)} latest seen")
            
            index1 = all_rc_graphs[rc1.perm][len(rc1)][rc1]
            index2 = all_rc_graphs[rc2.perm][len(rc2)][rc2]
            print(f"{index1} * {index2} = {coeff}")
            try:
                assert (latest_seen_left[perm][len(rc1)].get(rc1.perm,-1) < index1) or (latest_seen_right[perm][len(rc2)].get(rc2.perm,-1) < index2)
            except AssertionError:
                print(f"Order violation: {perm.trimcode}: {rc1.perm.trimcode} {len(rc2.perm.trimcode)}")
                print("********VIOLATION*******")
                print(f"  {latest_seen_left[perm][len(rc1)].get(rc1.perm)} * {latest_seen_right[perm][len(rc2)].get(rc2.perm)} latest seen")
                print(f"  {index1} * {index2}")
                print("********END VIOLATION*******")
            #print(f"  {index1} * {index2} = {coeff}")
            # latest_seen_left[perm][len(rc1)][rc1.perm] = index1
            # latest_seen_right[perm][len(rc2)][rc2.perm] = index2
            for co_perm in co_principal_perms(rc1.perm):
                #print(f"    Co-principal: {co_perm.trimcode} {all_rc_graphs[co_perm][len(rc1)].get(rc1)}")
            #print(f"  {index1} * {index2} = {coeff}")
                if len(co_perm) <= n:
                    assert co_perm.trimcode < perm.trimcode
                    latest_seen_left[co_perm][len(rc1)][rc1.perm] = index1 if latest_seen_left[co_perm][len(rc1)].get(rc1.perm, -1) < index1 else latest_seen_left[co_perm][len(rc1)][rc1.perm]
                    latest_seen_right[co_perm][len(rc2)][rc2.perm] = index2 if latest_seen_right[co_perm][len(rc2)].get(rc2.perm, -1) < index2 else latest_seen_right[co_perm][len(rc2)][rc2.perm]
        print("---------")

if __name__ == "__main__":
    main()
