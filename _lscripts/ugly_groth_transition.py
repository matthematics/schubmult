from schubmult import *

r = RCGraphRing()
rw = WCGraphRing()

def groth_to_schub_as_rc(groth_perm: Permutation, length):
    ret = r.zero
    for perm, coeff in WCGraph.groth_to_schub(groth_perm, Gx._beta).items():
        ret += coeff*r.schub(perm, length)
    return ret

def schub_to_groth_as_rc(schub_perm: Permutation, length):
    ret = rw.zero
    for perm, coeff in WCGraph.schub_to_groth(schub_perm, Gx._beta).items():
        ret += coeff*rw.groth(perm, length)
    return ret

def main(n):
    perms = Permutation.all_permutations(n)
    for perm in perms:
        # groth_poly = Gx(perm).expand()
        # groth_test = rw.groth(perm).polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True).expand()
        # assert (groth_poly - groth_test).expand() == 0, f"Groth polynomial test failed for {perm} in S_{n}\n{groth_poly} != {groth_test}"
        length = max(1,len(perm.trimcode))
        schub1 = groth_to_schub_as_rc(perm, length)
        trans_schub = schub1.zero_out_last_row()
        groth_trans = rw.zero
        for rc, coeff in trans_schub.items():
            if rc.is_principal:
                groth_trans += coeff*schub_to_groth_as_rc(rc.perm, length - 1)
        test_groth = Gx.from_expr(Gx(perm).expand().subs(Sx.genset[length],0))
        for wc, coeff in groth_trans.items():
            diff = test_groth.get(wc.perm, 0) - coeff
            print(coeff)
            assert diff.expand() == 0, f"Groth polynomial test failed for {perm} in S_{n}\n{diff}"
        #schub_test = schub_to_groth_as_rc(perm, len(perm.trimcode)).polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True).expand()
        #assert (groth_trans - test_groth).expand() == 0, f"Groth polynomial test failed for {perm} in S_{n}\n{groth_trans} != {test_groth}"

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    main(n)