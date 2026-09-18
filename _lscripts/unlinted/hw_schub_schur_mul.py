from schubmult import *
from schubmult.utils.tuple_utils import pad_tuple

g = BoundedRCFactorAlgebra()
r = RCGraphRing()

def cem_schub(perm, n):
    return sum([g.from_tensor_dict(cem_dict, n) for rc, cem_dict in RCGraph.full_CEM(perm, n).items()])

def cem_schub_schur_decomp(perm, n):
    from sympy import Mul, Add, expand
    result = Sx.zero @ Sx.zero
    cd = perm.strict_mul_dominant().trimcode
    if any(a < n for a in cd):
        toadd = min(n - a for a in cd if a < n)
        cd = [a + toadd for a in cd]
    domperm = uncode(cd)
    reppy = expand(Sx(perm).cem_rep(mumu=~domperm, elem_func=Sx.symbol_elem_func), func=False)
    for arg in Add.make_args(reppy):
        coeff, schur_part = arg.as_coeff_Mul()
        part1 = Sx.one
        part2 = Sx.one
        for elem_arg in Mul.make_args(schur_part):
            if elem_arg.numvars < n:
                part1 *= elem_arg
            else:
                part2 *= elem_arg
        result += coeff * part1 @ part2
    return result

if __name__ == "__main__":
    from sympy import pretty_print
    import sys
    import itertools
    n = int(sys.argv[1])
    extra = 1
    perms = [perm for perm in Permutation.all_permutations(n + extra) if len(perm.trimcode) <= n]
    
    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        cem_elem1 = cem_schub_schur_decomp(perm1, n)
        cem_elem2 = cem_schub_schur_decomp(perm2, n)

        result = r.zero

        prd = Sx(perm1) * Sx(perm2)

        for (base_perm, grass_perm), coeff1 in cem_elem1.items():
            for (base_perm2, grass_perm2), coeff2 in cem_elem2.items():
                graph_base = (cem_schub(base_perm, n)*cem_schub(base_perm2, n)).to_rc_graph_ring_element().resize(n)
                graph_grass = ((cem_schub(grass_perm, n) *cem_schub(grass_perm2, n))).to_rc_graph_ring_element().resize(n)
                for rc_grass, coeff4 in graph_grass.items():
                    result_base = r.zero
                    for rc, coeff3 in graph_base.items():                    
                            if rc.is_principal:
                                assert coeff3 >= 0, f"Negative coefficient for {rc} in product of {perm1} and {perm2} with base perms {base_perm}, {base_perm2} and grass perms {grass_perm}, {grass_perm2}"
                                result_base += coeff3 * r(rc.to_highest_weight()[0])
                    for rc_base, coeff3 in graph_base.items():
                        rcc = rc_base.squash_product(rc_grass)
                        if rcc.is_highest_weight and rcc.extremal_weight == pad_tuple(rcc.perm.trimcode, len(rcc)): 
                            result += coeff1 * coeff2 * coeff3 * coeff4 * r(rcc)

        prd2 = Sx.zero
        for rc, coeff in result.items():
            assert coeff == prd.get(rc.perm, 0), f"Failure for {perm1}, {perm2} at {rc.perm} with coeff {coeff} vs {prd.get(rc.perm, 0)}"
            prd2 += coeff * Sx(rc.perm)

        assert prd == prd2, f"Failure for {perm1}, {perm2} with result {result}, expected {prd}, got {prd2}"
        print("Success for", perm1, perm2)

        
        
