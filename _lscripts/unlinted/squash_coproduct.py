from schubmult import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    rcs = []
    grass_rcs = []
    for i in range(n):
        rcs.append(set())
        grass_rcs.append(set())

    for perm in perms:
        k = len(perm.trimcode)
        rc_set = RCGraph.all_rc_graphs(perm, k)
        rcs[k].update(rc_set)
        if len(perm.descents()) == 1:
            grass_rcs[k].update(rc_set)
    coprod_dict = {}
    R = RCGraphRing()
    tring = R @ R
    for max_d in range(n):
        for g_rc in grass_rcs[max_d]:
            for j in range(max_d + 1):
                for rc in rcs[j]:
                    rc2 = rc.resize(max_d)
                    the_rc = rc2.squash_product(g_rc)
                    coprod_dict[the_rc] = coprod_dict.get(the_rc, tring.zero) + tring((rc2, g_rc))
    print("Coproduct pairs (squash coproduct):")
    for rc, pickle in coprod_dict.items():
        for rc2, pickle2 in coprod_dict.items():
            the_stick = pickle * pickle2
            pigeon = rc * rc2
            for (suck, fish), coeff0 in the_stick.items():
                if len(fish.perm.trimcode) == len(fish):
                    found = False
                    pants = suck.squash_product(fish)
                    for suck2, coeff in pigeon.items():
                        if suck2 == pants:
                            found = True
                            break
                    if not found:
                        print(f"Missing pair: {suck} * {fish} = {pants} not in {pigeon}")
                        exit(1)
                    else:
                        print("MOSFO")
