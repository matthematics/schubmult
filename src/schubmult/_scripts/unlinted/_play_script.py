from schubmult import Permutation, elem_sym_perms, elem_sym_perms_q, q_vector, phi_d, tau_d, RCGraph, SchubertBasis, ASx, WordBasis, FA, Sx, RCGraphRing, Plactic, RootTableau

if __name__ == "__main__":
    import itertools
    n = 6
    perms = Permutation.all_permutations(n)

    mn = 10000

    for perm1, perm2 in itertools.combinations_with_replacement(perms, 2):
        if any(v > 1 for v in (Sx(perm1) * Sx(perm2)).values()) and max(len(perm1.descents()), len(perm2.descents())) > 1:
            if max(len(perm1.trimcode),len(perm2.trimcode))*max(perm1.inv, perm2.inv) <= mn:
                print(f"{perm1} {perm2} inv={max(len(perm1.trimcode),len(perm2.trimcode))*max(perm1.inv, perm2.inv)} Sx1*Sx2={Sx(perm1) * Sx(perm2)}")
                mn = max(len(perm1.trimcode),len(perm2.trimcode))*max(perm1.inv, perm2.inv)