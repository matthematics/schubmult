from schubmult import *
from sympy import pretty_print
from schubmult.rings.free_algebra import *

EE = FreeAlgebra(ElementaryBasis)

def rc_to_ee(rc):
    if isinstance(rc, RCGraph):
        # rc_inv = rc.transpose(len(rc.perm))
        # vec = list(reversed(rc_inv.length_vector))
        # vec = vec[:len(rc) - 1] + sorted(vec[len(rc) - 1:], reverse=True)
        
        # #vec = list(reversed([*rc_ni.length_vector]))
        # while len(vec)>0 and vec[-1] == 0:
        #     vec.pop()
        # result = EE(tuple(vec),len(rc))
        # rc_inv = rc_inv.resize(len(rc.perm) - 1)
        # vec2 = list(reversed(rc_inv.length_vector))
        # vec2 = vec2[:len(rc) - 1] + sorted(vec2[len(rc) - 1:], reverse=True)
        
        # #vec = list(reversed([*rc_ni.length_vector]))
        # while len(vec2)>0 and vec2[-1] == 0:
        #     vec2.pop()
        # vec = tuple(vec)
        # vec2 = tuple(vec2)
        # if vec == vec2:
        #     return result
        #return result + EE(tuple(vec2),len(rc))
        
        pants = ASx(rc.perm, len(rc)).change_basis(ElementaryBasis)
        pants2 = FA(*rc.length_vector).change_basis(ElementaryBasis)
        # def trim(vec):
        #     vec = list(vec)
        #     vec = vec[:len(rc) - 1] + sorted(vec[len(rc) - 1:], reverse=True)
        #     while len(vec)>0 and vec[-1] == 0:
        #         vec.pop()
        #     return tuple(vec[:len(rc) - 1] + sorted(vec[len(rc) - 1:], reverse=True))
        # return sum([coeff * EE(trim(key1), len(rc)) for (key1, key2), coeff in pants.items()])
        result = EE.from_dict({key: pants[key] for key, coeff in pants2.items() if key in pants})
        # pretty_print(result)
        #return FA(*rc.length_vector) @ result
        return result
    result = 0
    for rc0, coeff in rc.items():
        result += coeff * rc_to_ee(rc0)
    return result

r = RCGraphRing()

def perm_length_to_rc(perm, length_vector):
    return r.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, len(length_vector), weight=length_vector), 1))

def rc_to_pl(rc):
    if isinstance(rc, RCGraph):
        return ASx(rc.perm, len(rc))@FA(*rc.length_vector)
    return sum([coeff * rc_to_pl(rc0) for rc0, coeff in rc.items()])

def pare(tensor):
    result = 0
    for (perm_elem1, vec1, perm_elem2, vec2), coeff in tensor.items():
        if FA(*vec1).change_basis(SchubertBasis).get(perm_elem1, 0) != 0 and FA(*vec2).change_basis(SchubertBasis).get(perm_elem2, 0) != 0:
            result += coeff * tensor.ring((perm_elem1, vec1, perm_elem2, vec2))
    return result



if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])

    
    perms = Permutation.all_permutations(n)

    
    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     farpo = set()
    #     for rc in RCGraph.all_rc_graphs(perm, len(perm) - 1):
    #         to_add = tuple(sorted(tuple(rc_to_ee(rc).items())))
    #         if to_add not in farpo:
    #             #to_add = rc_to_ee(rc)
    #             farpo.add(to_add)
    #         else:
    #             print("Duplicate found for", perm.trimcode, "with RC graph") 
    #             pretty_print(rc)
                
            #spitoon = rc_to_ee(RCGraph.principal_rc(perm, len(perm.trimcode)))
        # if perm.inv == 0:
        #     continue
        # rc_prod = sum([rc_to_ee(rc) for rc in RCGraph.all_hw_rcs(perm, n - 1)]).change_basis(SchubertBasis)
        # assert rc_prod.almosteq(ASx(perm,n-1)), f"Failure for {perm.trimcode}, got {rc_prod} instead of {ASx(perm, len(perm.trimcode))}"
    # is a sub coalgebra
    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        the_prod = Sx(perm1) * Sx(perm2)
        length = max(len(perm1.trimcode), len(perm2.trimcode))
        result = 0
        prod = itertools.product(RCGraph.all_rc_graphs(perm1, length), RCGraph.all_rc_graphs(perm2, length))
        # key_set = set()
        # for rc1, rc2 in prod:
        #     if (rc1, rc2) in key_set or (rc2, rc1) in key_set:
        #         continue
        #     key_set.add((rc1, rc2))
        # dup_key_set = set(key_set)
        used = set()
        
            # rc3 = RCGraph.principal_rc(perm3, length)
        old_pairs = set(prod)
        for perm3, coeff in the_prod.items():
            pairs = set(old_pairs)
            
            
            
            # if combined_length != rc3.length_vector:
            #     continue
            
            for rc3 in RCGraph.all_rc_graphs(perm3, length):
                for rc1,rc2 in old_pairs:
                    combined_length = tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)])
                    if rc3.length_vector != combined_length:
                        continue
                    rc3_coprod = rc_to_pl(rc3).coproduct()        
                    # elem_prod = rc_to_pl(rc1) * rc_to_pl(rc2)
                    # rc_elem_prod = rc_to_pl(r(rc1)* r(rc2))
                    # assert elem_prod.almosteq(rc_elem_prod), f"Failure for {perm1.trimcode}, {perm2.trimcode}, got {elem_prod} instead of {rc_elem_prod}\n{rc1}\n{rc2}\n{r(rc1)*r(rc2)}"
                    # rc1_coprod = rc_to_pl(rc1).coproduct()
                    # rc2_coprod = rc_to_pl(rc2).coproduct()

                    # product_coprod = rc1_coprod * rc2_coprod

                    # rc_prod = r(rc1) * r(rc2)

                    # rc_product_coprod = rc_to_pl(rc_prod).coproduct()

                    # assert product_coprod.almosteq(rc_product_coprod), f"Failure for {perm1.trimcode}, {perm2.trimcode}, got {product_coprod} instead of {rc_product_coprod}\n{rc1_coprod}\n{rc2_coprod}\n{rc_prod}\n{rc_product_coprod-product_coprod}"
                    
                        #for rc3 in RCGraph.all_rc_graphs(perm3, length, weight=tuple([a + b for a, b in zip(rc1.length_vector, rc2.length_vector)])):
                    
                    
                    if rc3_coprod.get(((rc1.perm, length), rc1.length_vector, (rc2.perm, length), rc2.length_vector), 0) != 0:# and (rc1, rc2) not in used:
                        result +=  r(rc3)
                        pairs.remove((rc1, rc2))
                        #break
                old_pairs = set(pairs)
                        
            
                    
                #break
        #pretty_print(result)
        the_result = Sx.from_dict({rc.perm: coeff for rc, coeff in result.items() if rc.is_principal})
        #the_result = Sx.from_expr(result)#Sx.from_dict({k: v for (k, _), v in result.to_free_algebra_element().items()})
        assert the_result == the_prod, f"Failure for {perm1.trimcode}, {perm2.trimcode}, instead of {the_prod}\ngot\n{the_result}"
        print("Fated successfully for", perm1.trimcode, perm2.trimcode)
        #r(rc1)*r(rc2)
        # sumup = rc_to_ee(rc_prod)
        # assert prodo.almosteq(sumup), f"Failure for {perm1.trimcode}, {perm2.trimcode}, got {prodo} instead of {sumup}"