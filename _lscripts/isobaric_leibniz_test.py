from schubmult import *

def apply_isobaric(i, schub_elem, beta):
    return (1 + beta * schub_elem.ring.genset[i]) * schub_elem.divdiff(i) - beta * schub_elem


if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for k in range(1, n):
        for perm1, perm2 in itertools.product(perms, repeat=2):
            the_prod = Sx(perm1) * Sx(perm2)

            the_prod_for_real = apply_isobaric(k, the_prod, Gx._beta) + Gx._beta * the_prod

            leibniz_prod = (apply_isobaric(k, Sx(perm1), Gx._beta) +Gx._beta * Sx(perm1)) * Sx(perm2).simpleref(k)
            leibniz_prod += Sx(perm1) * (apply_isobaric(k, Sx(perm2), Gx._beta) + Gx._beta * Sx(perm2))

            assert leibniz_prod.almosteq(the_prod_for_real), f"Leibniz rule failed for {perm1}, {perm2}, {k}: \n{the_prod_for_real=}\n{leibniz_prod=}"

