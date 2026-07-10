from schubmult.rings.free_algebra import *
from schubmult import *
from sympy import sympify, pretty_print

if __name__ == "__main__":
    import sys



    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.pad_code(n - 1)) for perm in perms]
    GroveDual = FreeAlgebra(GroveBasis)
    for u in comps:
        
        for v in comps:
            # if v.inv == 0:
            #         continue
            for i in range(len(u), n):
            
                
                for j in range(len(v), n):
                    the_product = GroveDual(*pad_tuple(u, i)) * GroveDual(*pad_tuple(v, j))
                    # print(f"Product of {(u, i), (v, j)} = {the_product}")
                    #result_dict = {}
                    the_uc = next(iter([rc for rc in WCGraph.grove_wcs(u, i) if rc.is_reduced]))
                    if len(the_uc) != i:
                        raise ValueError(f"Unexpected length {len(the_uc)} for u={u} and i={i}")
                    
                    the_vc = next(iter([rc for rc in WCGraph.grove_wcs(v, j) if rc.is_reduced]))
                    if len(the_vc) != j:
                        raise ValueError(f"Unexpected length {len(the_vc)} for v={v} and j={j}")
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
                        
                        for wc in WCGraph.grove_wcs(w, i + j):
                            up, down = wc.vertical_cut(i)
                            if len(up) != i or len(down) != j:
                                raise ValueError(f"Unexpected lengths {len(up)}, {len(down)} for vertical cut of {wc} at {i}")
                            if up == the_uc and down == the_vc:
                                result += 1
                        try:
                            assert result == coeff, f"Mismatch for {(u, i), (v, j)}: expected {coeff} for {w=}, got {result}"
                        except AssertionError as e:
                            print(f"FAIL: {(u, i), (v, j)}: expected {coeff}, got {result}")
                            pretty_print(the_uc)
                            pretty_print(the_vc)
                            for w in WCGraph.all_wc_graphs(w, length):
                                pretty_print(w.vertical_cut(i))
                            raise e
                    print(f"I love paint {u=} {v=}")