from schubmult import *

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    for w in perms:
        fordict = {}
        for wc in WCGraph.all_wc_graphs(w, n - 1):
            word = []
            compat = []
            perm_last = Permutation([])
            letter_last = 0
            append_val = 1
            for a in wc.perm_word:
                new_perm = perm_last @ Permutation.ref_product(a)
                if new_perm != perm_last:
                    word.append(a)
                    perm_last = new_perm
                    if a >= letter_last:
                        append_val += 1
                    compat.append(append_val)
                    letter_last = a
            rc = RCGraph.from_reduced_compatible(word, compat)
            fordict[rc.forest_invariant] = fordict.get(rc.forest_invariant, 0) + wc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True)
        testo = {}
        for invar, val in fordict.items():
            cd = invar.forest.code
            if cd in testo:
                assert (testo[cd] - val).expand() == 0, f"Mismatch for invariant {invar} code {cd}: {testo[cd]} != {val}"
            else:
                testo[cd] = val
        print("happy")