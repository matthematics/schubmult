from schubmult import *
from sympy import pretty_print

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    for perm1, perm2 in itertools.combinations(perms, 2):
        # if len(perm1.descents()) > 1 or len(perm2.descents()) > 1 or max(perm1.code, default=0) > 1 or max(perm2.code, default=0) > 1:
        #     continue
        
        # if perm1.is_dominant or perm2.is_dominant:
        #     continu
    
        for rc1, rc2 in itertools.product(RCGraph.all_hw_rcs(perm1, n - 1), RCGraph.all_hw_rcs(perm2, n - 1)):
            try:
                    #print("Trying to compute product in RC graph ring...")
                    r(rc1) % r(rc2)
                    #print("Success!")
                    continue
            except NotImplementedError:
                pass
            weight1, dperm1 = rc1.classify_demazure_crystal()
            weight2, dperm2 = rc2.classify_demazure_crystal()
            # Permutation.does_demazure_crystal_tensor_decompose(weight1, dperm1, weight2, dperm2)
            
            
            tryit = Permutation.does_demazure_crystal_tensor_decompose(weight1, dperm1, weight2, dperm2)
            if not tryit:
                continue
            pretty_print(r(rc1) @ r(rc2))
            print(f"Demazure crystal of {perm1} with dominant weight {weight1} and lowest weight {dperm1}")
            print(f"Demazure crystal of {perm2} with dominant weight {weight2} and lowest weight {dperm2}")
            print()
            print(f"Decomposes? {tryit}")
            print(f"Other direction? {Permutation.does_demazure_crystal_tensor_decompose(weight2, dperm2, weight1, dperm1)}")
            if tryit:                
                    produc = Sx(rc1.perm) * Sx(rc2.perm)
                    rc1_t = rc1.transpose().normalize().to_highest_weight()[0]
                    rc2_t = rc2.transpose().normalize().to_highest_weight()[0]
                    combinat = r(rc2_t) * r(rc1_t)
                    gotone = False
                    for rc, coeff in combinat.items():
                        rc = rc.transpose(len(rc1))
                        
                        if produc.get(rc.perm, 0) != 0:# and rc.length_vector == tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)]):
                            gotone = True
                            pretty_print((rc, coeff))
                            print(f"Product of {rc1} and {rc2} has term {rc} with coefficient {coeff}, which corresponds to {Sx(rc.perm)} in the Schubert basis, which has product {produc.get(rc.perm, 0)}")
                    if not gotone:
                        print("WAHBOOGER")