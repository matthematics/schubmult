from schubmult import *
from schubmult.combinatorics.pipe_dream import PipeDream
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import *
from schubmult.combinatorics.indexed_forests import *

def groth_to_schub_as_rc(groth_perm: Permutation):
    boip = WCGraph.all_wc_graphs(groth_perm, len(groth_perm))
    #spits = [PipeDream.from_rc_graph(bpd).co_pipe_dream() for bpd in boip if PipeDream.from_rc_graph(bpd).co_pipe_dream().is_reduced]
    ret = {}
    for wc in boip:
        print(wc)
        cpdb = PipeDream.from_rc_graph(wc).co_pipe_dream()
        print(cpdb)
        print(cpdb.perm)
        print(cpdb.perm_word)
        print(cpdb.is_reduced)
        if cpdb.is_reduced:
            w0 = Permutation.w0(cpdb.rows)
            perm = cpdb.perm * w0
            
            ret[cpdb.co_pipe_dream()] = RCGraph.all_rc_graphs(perm)
    return ret


if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.pad_code(n - 1)) for perm in perms]
    # for k in range(1, n):
    #     for p in range(1, k + 1):
    for comp in comps:
        #    w = uncode([0] * (k - p) + [1] * p)
        dct = WCGraph.grove_to_forest(comp, beta=Gx._beta)
        # for bpd, rc_set in dct.items():
        #     assert len([bpd0 for bpd0 in BPD.all_unreduced_bpds(w) if bpd0.co_bpd().perm == bpd.co_bpd().perm]) == len(rc_set), f"Mismatch for {w}: {bpd=}, {len([bpd0 for bpd0 in BPD.all_unreduced_bpds(w) if bpd0.co_bpd().perm == bpd.co_bpd().perm])=}, {len(rc_set)=}"
        print(f"{comp=}: {dct}")
        poly1 = grove_polynomial(comp, Sx.genset, beta=Gx._beta)
        poly2 = 0
        for comp2, coeff in dct.items():
            poly2 += coeff * Forest(*comp2).expand()
        #poly2 = sum([rc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for bpd, rcs in groth_to_schub
        assert (poly1 - poly2).expand() == 0, f"Failed for {comp=}: {poly1=}, {poly2=}"
        print("Success ", comp)