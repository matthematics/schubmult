from schubmult import *
from symengine import Symbol
J = FreeAlgebra(JBasis)
t = Symbol("t")

def speed_pieri(p, j_elem):
    sch_elem = sum([coeff * ASx(uncode(c)) for c, coeff in j_elem.items()])
    up_elem = ASx(uncode([p])) * sch_elem 
    ret_elem = J.zero
    for (perm, n), coeff in up_elem.items():
        if 0 in perm.trimcode:
            ret_elem += t * coeff * J(*[a for a in perm.trimcode if a != 0])
        else:
            ret_elem += coeff * J(*perm.trimcode)
    return ret_elem

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

    for perm in perms:
        sputnik = perm.trimcode
        aplach = [a for a in sputnik if a != 0]
        jj = J(*aplach)
        for p in range(1, n):
            test1 = J(p) * jj
            test2 = speed_pieri(p, jj)
            assert test1 == test2, f"Failed for {perm} and p={p}: {test1} != {test2}"
    print("All tests passed!")