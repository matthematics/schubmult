from schubmult.rings.free_algebra import *
from schubmult.rings.free_algebra.grove_basis import GroveBasis
from schubmult import *
from schubmult.utils.tuple_utils import pad_tuple
from sympy import sympify, pretty_print

if __name__ == "__main__":
    import sys



    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.trimcode) for perm in perms]
    GroveDual = FreeAlgebra(GroveBasis)
    for u in comps:
        
        for v in comps:
            # if v.inv == 0:
            #         continue
            for i in range(len(u), n):
            
                padded_u = pad_tuple(u, i)        
                for j in range(len(v), n):
                    padded_v = pad_tuple(v, j)
                    the_product = GroveDual(*padded_u) * GroveDual(*padded_v)
                    # print(f"Product of {(u, i), (v, j)} = {the_product}")
                    #result_dict = {}
                    the_uc = set()
                    the_vc = set()
                    for (perm, _), cofo in GroveDual(*padded_u).change_basis(GrothendieckBasis).items():
                        if cofo == 0:
                            continue
                        the_uc.update({wc.grove_invariant for wc in WCGraph.all_wc_graphs(perm, i, weight=padded_u) if wc.grove_weight == padded_u})

                    for (perm, _), cofo in GroveDual(*padded_v).change_basis(GrothendieckBasis).items():
                        if cofo == 0:
                            continue
                        the_vc.update({wc.grove_invariant for wc in WCGraph.all_wc_graphs(perm, j, weight=padded_v) if wc.grove_weight == padded_v})
                                
                    
                    for w, coeff in the_product.items():
                        # print(f"Checking {(w, length)} with coeff {coeff}")
                        if len(w) > i + j:
                            print("waring")
                            continue
                        
                        result = 0
                        #coeff = sympify(coeff).subs(Gx._beta, -1)
                        if coeff == 0:
                            continue
                        # if length != i + j:
                        #     print(f"Unexpected length {length} for product of {(u, i), (v, j)}, {coeff=}")
                        #     continue
                        gwdct = {}
                        ginv = None
                        for wc in WCGraph.all_wc_graphs(uncode(w), i + j):
                            if wc.grove_weight == w:
                                ginv = wc.grove_invariant
                                break
                        
                        for wc in WCGraph.grove_wcs(w, i + j, ginv):
                            up, down = wc.vertical_cut(i)
                            if len(up) != i or len(down) != j:
                                raise ValueError(f"Unexpected lengths {len(up)}, {len(down)} for vertical cut of {wc} at {i}")
                            if up in the_uc and down in the_vc:
                                result += 1
                        try:
                            assert result == coeff, f"Mismatch for {(padded_u, i), (padded_v, j)}: expected {coeff} for {w=}, got {result}"
                        except AssertionError as e:
                            print(f"FAIL: {(padded_u, i), (padded_v, j)}: {w=} expected {coeff}, got {result}")
                            print("WNATU")
                            pretty_print(the_uc)
                            print("WNATV")
                            pretty_print(the_vc)
                            print("BUGS")
                            for w_ in WCGraph.all_wc_graphs(uncode(w), i + j):
                                pretty_print(w_.vertical_cut(i))
                            raise e
                    print(f"I love paint {padded_u=} {padded_v=}")