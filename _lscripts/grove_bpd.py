from schubmult import *
from schubmult.rings.polynomial_algebra import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.symbolic.poly.schub_poly import *
from schubmult.symbolic.poly.variables import *

def groth_to_schub_as_rc(groth_perm: Permutation):
    boip = BPD.all_unreduced_bpds(groth_perm, len(groth_perm))
    spits = [bpd.co_bpd() for bpd in boip if bpd.co_bpd().is_reduced]
    ret = {}
    for cpdb in spits:
        perm = cpdb.perm * Permutation.w0(cpdb.rows)
        ret[cpdb.co_bpd()] = RCGraph.all_rc_graphs(perm)
    return ret

def schub_to_groth(the_perm, beta=Symbol("beta")):
    for bpd in RCGraph.all_rc_graphs(the_perm, len(the_perm)):
        pd = PipeDream.from_rc_graph(bpd).co_pipe_dream()#.perm * w0
        rc = pd.to_wc_graph()
        dct[comp] = dct.get(bpd, S.Zero) + (-beta) ** (permo.inv - the_perm.inv)

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    zz = ZeroGeneratingSet()
    perms = Permutation.all_permutations(n)
    CompSchub = PolynomialAlgebra(CompositionSchubertPolyBasis(Sx.genset))
    for w in perms:
            #dct = schub_to_groth(w)
            # for bpd, rc_set in dct.items():
            #     assert len([bpd0 for bpd0 in BPD.all_unreduced_bpds(w) if bpd0.co_bpd().perm == bpd.co_bpd().perm]) == len(rc_set), f"Mismatch for {w}: {bpd=}, {len([bpd0 for bpd0 in BPD.all_unreduced_bpds(w) if bpd0.co_bpd().perm == bpd.co_bpd().perm])=}, {len(rc_set)=}"
            # print(f"{w.trimcode}: {dct}")
            comp = tuple(w.trimcode)
            poly1 = grove_polynomial(comp, Sx.genset, Gx._beta)#Forest(*comp).expand()#, Sx.genset, Gx._beta)
            forp = Forest.from_expr(poly1, length=len(comp))
            poly2 = grothendieck_poly(w, Sx.genset, zz, Gx._beta)
            forp2 = CompSchub.from_expr(poly2, length=len(comp))
            forp3 = Forest.from_expr(poly2, length=len(comp))
            print(f"Grove {comp}: {forp}")
            print(f"Groth {comp}: {forp2}")
            print(f"Grovh {comp}: {forp3}")
            # for bpd, val in dct.items():
            #     poly2 += val * grove_polynomial(tuple(rc.forest_weight), Sx.genset, Gx._beta)
            # #poly2 = sum([rc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for bpd, rcs in groth_to_schub
            # assert (poly1 - poly2).expand() == 0, f"Failed for {w}: {poly1=}, {poly2=}"
            # print("Success ", w)