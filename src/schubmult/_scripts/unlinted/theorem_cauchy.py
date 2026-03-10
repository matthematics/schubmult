from schubmult import *
from schubmult.utils.perm_utils import artin_sequences

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    r = RCGraphRing()
    dr = DualRCGraphRing()
    first_elem = 0
    second_elem = 0
    # for perm in perms:
    #     first_elem += dr.schub(~perm, n) @ dr.schub(perm, n)
    seqs = artin_sequences(n - 1)
    for seq in seqs:
        term = r.one @ dr.one.resize(n - 1)
        #assert len(seq) == n - 1
        for index, i in enumerate(tuple(reversed(seq)), start=1):
            #second_elem += r.sx(i, n) @ dr.sx(i, n)
            # panties = dr.full_elem_sym(i, index, n)
            # for rc in panties:
            #term = term * (Sx.from_expr(Sx.genset[index]**i) @ (RCGraph.one_row(i) ) @ dr.full_elem_sym(i, index, n))
            #term = term * (Sx.from_expr(Sx.genset[n - index]**(index - i)) @ dr.full_elem_sym(i, index, n))
            term = (r(RCGraph.one_row(i)) @ dr.one.resize(n - 1)) * term * (r.one @ dr.full_elem_sym(i, index, n - 1))
        second_elem += term
    second_elem = second_elem.ring.from_dict({(rc1, rc2): coeff for (rc1, rc2), coeff in second_elem.items() if coeff != 0 and len(rc1.perm) <= n})# and perm == rc1.perm})
    #assert first_elem.almosteq(second_elem), f"Failure: first_elem={first_elem}, second_elem={second_elem}"
    print(second_elem)
        