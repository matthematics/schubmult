from schubmult import *

if __name__ == "__main__":
    from schubmult.symbolic import expand, S
    from schubmult.abc import E
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    x = Sx.genset
    for perm in perms:
        if perm.inv == 0:
            continue
        dual_schub = ASx(perm, n - 1).change_basis(WordBasis)
        sm = 0
        mx = -5
        for word in dual_schub:
            m = max([a + i for i, a in enumerate(word)])
            if m > mx:
                mx = m
        #mx = n
        for word, coeff in dual_schub.items():
            elem_prod = S.One
            for i, a in enumerate(word):
                elem_prod *= E(mx - i - a, mx - i, x[1:], [0 for _ in range(3*n)])
            sm += coeff * elem_prod
        schub_elem = Sx.from_expr(sm)
        print(f"{perm}: {dual_schub} => {schub_elem}")
        #print(schub_elem)
        # assert all(coeff >= 0 for coeff in schub_elem.values()), f"Failed for {perm}, got {schub_elem}"
        # assert len({k: v for k, v in schub_elem.items() if v != 0}) == 1, f"Failed for {perm}, got {schub_elem}"