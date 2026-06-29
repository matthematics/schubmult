from schubmult import *

def groth_to_schub_as_rc(groth_perm: Permutation):
    boip = BPD.all_unreduced_bpds(groth_perm, len(groth_perm))
    spits = [bpd.co_bpd() for bpd in boip if bpd.co_bpd().is_reduced]
    ret = {}
    for cpdb in spits:
        perm = cpdb.perm * Permutation.w0(cpdb.rows)
        ret[cpdb.co_bpd()] = RCGraph.all_rc_graphs(perm)
    return ret


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    for k in range(1, n):
        for p in range(1, k + 1):
            w = uncode([0] * (k - p) + [1] * p)
            dct = groth_to_schub_as_rc(w)
            # for bpd, rc_set in dct.items():
            #     assert len([bpd0 for bpd0 in BPD.all_unreduced_bpds(w) if bpd0.co_bpd().perm == bpd.co_bpd().perm]) == len(rc_set), f"Mismatch for {w}: {bpd=}, {len([bpd0 for bpd0 in BPD.all_unreduced_bpds(w) if bpd0.co_bpd().perm == bpd.co_bpd().perm])=}, {len(rc_set)=}"
            print(f"{w.trimcode}: {dct}")
            poly1 = Gx(w).as_polynomial()
            poly2 = 0
            for bpd, rcs in dct.items():
                poly2 += sum([(Gx._beta ** (rc.perm.inv - w.inv)) * rc.polyvalue(Sx.genset) for rc in rcs])
            #poly2 = sum([rc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for bpd, rcs in groth_to_schub
            assert (poly1 - poly2).expand() == 0, f"Failed for {w}: {poly1=}, {poly2=}"
            print("Success ", w)