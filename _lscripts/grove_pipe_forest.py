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
    # comps = [tuple(perm.pad_code(n - 1)) for perm in perms]
    # # for k in range(1, n):
    # #     for p in range(1, k + 1):
    Grove1 = PolynomialAlgebra(GrovePolyBasis(Gx.genset, 1))
    for perm in perms:
        #    w = uncode([0] * (k - p) + [1] * p)
        dct = {}
        for wc in WCGraph.all_wc_graphs(perm, n - 1):
            wt = wc.grove_invariant
            dct[wt] = dct.get(wt, set())
            dct[wt].add(wc)
        for wt, wcs in dct.items():
            the_sum = sum([wc.polyvalue(Sx.genset, beta=1) for wc in wcs])
            grove_poly = Grove1(*wt.length_vector).expand()
            assert (the_sum - grove_poly).expand() == 0, f"Mismatch for {perm}, {wt}: {the_sum} != {grove_poly}"
            print("GOOD", perm, wt)