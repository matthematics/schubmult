from schubmult import *

def transition_rc(rc):
    ring = RCGraphRing()
    if rc.perm.inv == 0:
        return (rc,)
    if len(rc[-1]) == 0:
        return (*transition_rc(rc.zero_out_last_row()), rc)
    lower_rc = rc.extend(1).exchange_property(len(rc)).zero_out_last_row()
    return (*transition_rc(lower_rc), rc) 


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for perm in perms:
        rc_set = RCGraph.all_rc_graphs(perm, len(perm.trimcode))
        for rc in rc_set:
            print(f"RC Graph for permutation {list(perm)}:")
            print(rc)
            print("Transition:")
            for tr in transition_rc(rc):
                print(tr)
            print()