from schubmult import *

if __name__ == "__main__":
    import sys
    import itertools
    n = int(sys.argv[1])
    g = BoundedRCFactorAlgebra()
    bad_patterns = [[4,1,3,2],[1,4,3,2],[3,1,4,2]]
    avoids_bad = lambda p: all(not p.has_pattern(pat) for pat in bad_patterns)
    #bad_patterns = [[3,2,1],[3,1,4,2]]

    perms = Permutation.all_permutations(n)
    
    # for perm in perms:
    #     if avoids_bad(perm):
    #         stack = [perm]
    #         while stack:
    #             the_perm = stack.pop()
    #             if not avoids_bad(the_perm):
    #                 #print(f"Found positive pattern-avoidance counterexample: {the_perm}")
    #                 raise ValueError(f"Found positive pattern-avoidance counterexample: {perm=} {the_perm=}")
    #             for d in (~the_perm).descents():
    #                 stack.append(~((~the_perm).swap(d, d+1)))
    r = RCGraphRing()
    #for perm1, perm2 in itertools.product(perms, repeat=2):
    for perm in perms:
        #if avoids_bad(perm1) and avoids_bad(perm2):
        # prd = Sx(perm1) * Sx(perm2)
        if avoids_bad(perm):
            # assert avoids_bad(perm1) and avoids_bad(perm2), f"Found positive pattern-avoidance counterexample: {perm1=} {perm2=} in product {prd}"
            # for perm in prd:
            for rc in RCGraph.all_rc_graphs(perm):
                # for j in range(len(rc3)):
                #     #if not all(avoids_bad(rc4.perm) for rc4 in  rc3.vertical_cut(j)):
                #     if not avoids_bad(rc3.vertical_cut(j)[1].perm):
                #         raise ValueError(f"Found positive pattern-avoidance counterexample: {perm1=} {perm2=} {perm=} {rc3.perm=} {j=}\n{rc3=}\n{rc3.vertical_cut(j)=}")
                copro0 = r(rc).coproduct()
                copro1_left = (r@r@r).zero
                for (rc1, rc2), coeff in copro0.items():
                    copro1_left += coeff * (r(rc1).coproduct() @ r(rc2))
                copro1_right = (r@r@r).zero
                for (rc1, rc2), coeff in copro0.items():
                    copro1_right += coeff * (r(rc1) @ r(rc2).coproduct())
                if not copro1_left.almosteq(copro1_right):
                    print(f"Found non-coassociative coproduct: {perm=} {rc=}\ncoproduct={copro0}\n(r(rc1).coproduct() @ r(rc2))={copro1_left}\n(r(rc1) @ r(rc2).coproduct())={copro1_right}\n{copro1_left - copro1_right=}")
                    sys.exit(1)
                # assert all([avoids_bad(rc1.perm),avoids_bad(rc2.perm)] for (rc1,rc2) in copro), f"Found positive pattern-avoidance counterexample: {perm1=} {perm2=} {perm=} {rc3.perm=}\n{rc3=}\ncoproduct={copro}"
                # assert all(v > 0 for v in copro.values()), f"Found non-positive coproduct coefficient: {perm1=} {perm2=} {perm=} {rc3.perm=}\n{rc3=}\ncoproduct={copro}"
    #print("All checks passed!")
            # for length1, length2 in itertools.product(range(len(perm1.trimcode), n), range(len(perm2.trimcode), n)):
            #     for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, length1), RCGraph.all_rc_graphs(perm2, length2)):
            #         rc_elem = r(rc1) * r(rc2)
            #         for rc3, coeff in rc_elem.items():
            #             if coeff != 0 and not avoids_bad(rc3.perm):
            #                 raise ValueError(f"Found positive pattern-avoidance counterexample: {perm1=} {perm2=} {rc3.perm=} coeff={coeff}")
            # for perm3, coeff in prd.items():
            #     if coeff != 0 and not avoids_bad(perm3):
            #         raise ValueError(f"Found positive pattern-avoidance counterexample: {perm1=} {perm2=} {perm3=} coeff={coeff}")
        #m = len(perm.trimcode) + 1
        # m = n
        # for length in range(len(perm.trimcode), m):
        #     schub = g.schub_elem(perm, len(perm.trimcode))
        #     schub = g.from_tensor_dict(schub, size=length)
        #     if any(perm.has_pattern(pat) for pat in bad_patterns):
        #             assert any(coeff < 0 for coeff in schub.values()), f"Failed to find negative coefficient for {perm} in length {length}"
        #     else:
        #         assert all(coeff >= 0 for coeff in schub.values()), f"Unexpected negative coefficient for {perm} in length {length}"
    print("All checks passed!")