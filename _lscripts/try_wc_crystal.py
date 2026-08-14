from schubmult import *

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
        seen = set()
        for rc in RCGraph.all_rc_graphs(perm):
            key = next(iter(bw.from_brc_elem(br.from_rc_graph(rc, n))))
            for cryselem in set([bw.make_key(key2, n) for key2 in key.full_crystal_bothways()]).difference(seen):
                seen.add(cryselem)
                if bw.key_to_wc_graph(cryselem).perm == perm:
                    if bw.key_to_wc_graph(cryselem).resize(perm.max_descent) in cachit:
                        print(f"Duplicate \n{cryselem=} for {perm=} in addition to\n{cachit[bw.key_to_wc_graph(cryselem).resize(perm.max_descent)]=}")                        
                    cachit[bw.key_to_wc_graph(cryselem).resize(perm.max_descent)] = cryselem
        for wc in WCGraph.all_wc_graphs(perm):
            if wc not in cachit:
                print(f"Missing {wc=}")
        print(len(cachit), "crystal elements for", perm, " vs ", len(WCGraph.all_wc_graphs(perm)), "WC graphs")
        
        #for wc in WCGraph.all_wc_graphs(perm):
            