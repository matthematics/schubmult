from schubmult import *

R = RCGraphRing()

def fa_to_rc(fa):
    result = R.zero
    for w, coeff in fa.items():
        res = R.one
        for a in w:
            res = res * R(RCGraph.one_row(a))
        result += coeff * res        
    return result

def rc_ext_mult(fa_tens):
    T = R @ R
    result = T.zero
    FA = FreeAlgebra(WordBasis)
    for (w1, w2), coeff in fa_tens.items():
        result += coeff * T.ext_multiply(fa_to_rc(FA(*w1)), fa_to_rc(FA(*w2)))
    return result

if __name__ == "__main__":
    import sys
    import itertools

    n = int(sys.argv[1])
    #for perm in Permutation.all_permutations(n):
        # for length in range(len(perm.trimcode), n):
        #     A = {}
        #     for rc in RCGraph.all_rc_graphs(perm, length):
        #         A[tuple(rc.length_vector)] = rc.polyvalue(Sx.genset)
        #     print(f"{perm} length {length}:")
        #     for w, coeff in A.items():
        #         print(f"  {w}: {coeff}")
        
        # A = ASx(perm, len(perm.trimcode)).change_basis(WordBasis).coproduct()
        # result = rc_ext_mult(A)

    for p in range(1, n + 1):
        elem = R.coproduct_on_basis(RCGraph.one_row(p))
        zz = R.coproduct_on_basis(RCGraph([()]))
        val = elem
        for _ in range(p):
            val *= zz
        print(val)
        for (rc1, rc2), coeff in val.items():
            assert (Sx(rc1.perm) * Sx(rc2.perm))