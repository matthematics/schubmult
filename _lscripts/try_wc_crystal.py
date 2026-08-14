from schubmult import *
from schubmult.rings.polynomial_algebra import *

br = BoundedRCFactorAlgebra()
bw = BoundedWCFactorAlgebra()

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    factorizations = {}
    
    for perm in perms:
        if perm.inv == 0:
            continue
        cachit = {}
        crystals = {}
        seen = set()
        # for rc in RCGraph.all_rc_graphs(perm):
        #     key = next(iter(bw.from_brc_elem(br.from_rc_graph(rc, n))))
        #     for cryselem in set([bw.make_key(key2, n) for key2 in key.full_crystal_bothways()]).difference(seen):
        #         seen.add(cryselem)
        #         if bw.key_to_wc_graph(cryselem).perm == perm:
        #             if bw.key_to_wc_graph(cryselem).resize(perm.max_descent) in cachit:
        #                 print(f"Duplicate \n{cryselem=} for {perm=} in addition to\n{cachit[bw.key_to_wc_graph(cryselem).resize(perm.max_descent)]=}")                        
        #             cachit[bw.key_to_wc_graph(cryselem).resize(perm.max_descent)] = cryselem
        # extra = {}
        signature = [(rc.perm.inv, len(rc)) for rc in next(iter(br.from_rc_graph(RCGraph.principal_rc(perm), n)))]
        all_cryselems = {bw.make_key(RCGraph([]), n)}
        for index in range(len(signature)):
            new_all_cryselems = set()
            for cryselem in all_cryselems:
                for wc in WCGraph.elem_sym_wcs(signature[index][0], signature[index][1]):
                    new_all_cryselems.add(bw.make_key((*cryselem,wc,), n))
            all_cryselems = new_all_cryselems
        for key in all_cryselems:
            if key in seen:
                continue
            #wc = bw.key_to_wc_graph(cryselem).resize(perm.max_descent)
            #crystal = set()
            for cryselem in set([bw.make_key(key2, n) for key2 in key.full_crystal_bothways(lambda x: bw.key_to_wc_graph(bw.make_key(x, n)).perm == perm)]).difference(seen):
                seen.add(cryselem)
                wc = bw.key_to_wc_graph(cryselem).resize(perm.max_descent)
                weights_seen = set()
                if wc.perm == perm:
                    # if bw.key_to_wc_graph(cryselem).resize(perm.max_descent) in cachit:
                    #     print(f"Duplicate \n{cryselem=} for {perm=} in addition to\n{cachit[bw.key_to_wc_graph(cryselem).resize(perm.max_descent)]=}")                        
                    if len(weights_seen) > 0:
                        if wc.extremal_weight not in weights_seen:
                            print(f"WAAAAAAH {perm}")
                    else:
                        weights_seen.add(wc.extremal_weight)
                    crystals.setdefault(wc.extremal_weight, set()).add(wc)
                    
                    #cachit[bw.key_to_wc_graph(cryselem).resize(perm.max_descent)] = cryselem
            # if len(crystal) > 0:
            #     crystals.append(crystal)
            
        # for wc in WCGraph.all_wc_graphs(perm):
        #     if wc not in cachit:
        #         print(f"Missing {wc=}")
        #         print("In extra?")
        #         if wc in extra:
        #             print("Yes, extra:", extra[wc])
        #         else:
        #             print("No, not in extra")
        #print(len(cachit)+len(extra), "crystal elements for", perm, " vs ", len(WCGraph.all_wc_graphs(perm)), "WC graphs")
        for weight, crystal in crystals.items():
            the_poly = sum([wc.polyvalue(Sx.genset) for wc in crystal])
            wc = next(iter(crystal))
            the_lascoux = LascouxPoly(*weight).expand()
            assert (the_poly - the_lascoux).expand() == 0, f"Error: Crystal sum for {perm} does not match Lascoux polynomial, got {(the_poly-the_lascoux).expand()} instead of 0"
        print("Painted bagels")
        
        #for wc in WCGraph.all_wc_graphs(perm):
            