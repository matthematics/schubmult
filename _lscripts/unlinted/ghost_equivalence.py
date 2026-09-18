from multiprocessing import Pool, cpu_count
from schubmult.rings.polynomial_algebra import *
from schubmult import *
from schubmult.abc import x
from schubmult.utils.perm_utils import mu_A
from schubmult.rings.free_algebra import *
from schubmult.utils.tuple_utils import pad_tuple
from schubmult.combinatorics.permutation import bad_classical_patterns
import itertools

def check_perm(perm_arr):
    perm = Permutation(perm_arr)
    n = len(perm_arr)
    r = RCGraphRing()
    for i in range(len(perm.trimcode), n):
        rc_elem = r.from_free_algebra_element(ASx(perm, i))
        if not all(int(v) in {1, 0, -1} for v in rc_elem.values()):
            return (False, perm_arr, dict(rc_elem))
    return (True, perm_arr, None)

EP = PolynomialAlgebra(ElemSymPolyBasis(Sx.genset))

# def in_good_algebra(elem_comp, monomial_vec):
#     #monom = EP.from_dict({elem_comp: 1}).change_basis(EP._basis.monomial_basis)
#     #return monom.get(monomial_vec, 0) != 0
#     assert isinstance(elem_comp[1], Permutation)
#     if FA(*)
#     return True

SS = FreeAlgebra(SchubertSchurBasis)
EE = FreeAlgebra(ElementaryBasis)


neg_pats = [Permutation([1, 4, 3, 2]), Permutation([4, 1, 3, 2]), Permutation([3, 1, 4, 2])]


g = GrassTensorAlgebra()
r = RCGraphRing()
def grass_tensor_elem(perm, n = None):
    if n is None:
        n = len(perm.trimcode)
    grass_elem = ring.zero
    cem = RCGraph.full_CEM(perm, n)
    
    for rc, cem_dict in cem.items():
        #if len(rc.perm) <= n:
        if True:
            grass_elem += g.from_dict(cem_dict)
        elif rc.perm != perm:
            grass_elem += g.from_dict(cem_dict)
        else:
            rc_base, rc_grass = rc.resize(n).squash_decomp()
            if rc_grass.perm.inv != 0:
                grass_elem += g.from_dict({(*k, rc_grass): v for k, v in RCGraph.full_CEM(rc_base.perm, n - 1)[rc_base.resize(n - 1)].items()})
            else:
                grass_elem += g.from_dict(cem_dict)
    return grass_elem

def n_tensor(key, n):
    factor1 = g.one
    factor2 = g.one
    flag = False
    for index in range(len(key)):
        if len(key[index].perm) <= n:
            factor1 *= g((key[index],))
        else:
            factor2 *= g((key[index],))

    return factor1 @ factor2

def doit(grass_elem1, n):
    picnic1 = (g@r).zero
    for key, coeff in grass_elem1.items():
        to_add1 = g.one
        index = 0
        for rc in key:
            if len(rc.perm.trimcode) < n - 1:
                to_add1 *= g((rc,))
                index += 1
            else:
                break
        assert (next(iter([len(rc.perm) for rc in to_add1.to_rc_graph_ring_element()]))) <= n - 1
        picnic1 += coeff * to_add1 @ g(key[index:]).to_rc_graph_ring_element()
    return picnic1

if __name__ == "__main__":
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.abc import *
    from schubmult.rings.polynomial_algebra import *
    from sympy import pretty_print
    import sys
    import itertools
    n = int(sys.argv[1])
    extra = 2
    perms = [perm for perm in Permutation.all_permutations(n + extra) if len(perm.trimcode) <= n]
    

    for perm1, perm2 in itertools.product(perms, repeat=2):
        cem1 = RCGraph.full_CEM(perm1, n)
        cem2 = RCGraph.full_CEM(perm2, n)

        # ghost1 = {rc: cem1[rc] for rc in cem1 if rc.perm != perm1}
        # ghost2 = {rc: cem2[rc] for rc in cem2 if rc.perm != perm2}

        # non_ghost1 = {rc: cem1[rc] for rc in cem1 if rc.perm == perm1}
        # non_ghost2 = {rc: cem2[rc] for rc in cem2 if rc.perm == perm2}

        # nontrivial_ghost1 = [((g.from_dict(g_cem1) * g.from_dict(elem2)).to_rc_graph_ring_element())  for g_cem1, elem2 in itertools.product(ghost1.values(), non_ghost2.values()) if not (g.from_dict(g_cem1) * g.from_dict(elem2)).to_rc_graph_ring_element().almosteq(r.zero)]
        # nontrivial_ghost2 = [((g.from_dict(elem1) * g.from_dict(g_cem2)).to_rc_graph_ring_element())  for g_cem2, elem1 in itertools.product(ghost2.values(), non_ghost1.values()) if not (g.from_dict(elem1) * g.from_dict(g_cem2)).to_rc_graph_ring_element().almosteq(r.zero)]
        # ghost_positive = r.zero
        # ghost_negative = r.zero
        # for ghost_elem in nontrivial_ghost1 + nontrivial_ghost2:
        #     for rc, coeff in ghost_elem.items():
        #         if coeff > 0:
        #             ghost_positive += coeff * r(rc)
        #         else:
        #             ghost_negative += -coeff * r(rc)
        
        result = r.zero
        prd = Sx(perm1) * Sx(perm2)
        for rc1, cem_dict1 in cem1.items():
            for rc2, cem_dict2 in cem2.items():
                cem_tensor1 = doit(cem_dict1, n + 1)#sum([v * doit(k, n) for k, v in cem_dict1.items()])
                cem_tensor2 = doit(cem_dict2, n + 1)#sum([v * doit(k, n) for k, v in cem_dict2.items()])
                # if (rc1.perm != perm1 or rc2.perm != perm2) and not (rc1.perm == perm1 and rc2.perm == perm2):
                #     ghost_elem = (g.from_dict(cem_dict1) * g.from_dict(cem_dict2)).to_rc_graph_ring_element()
                #     for rc, coeff in ghost_elem.items():
                #         if coeff > 0:
                #             ghost_positive -= coeff * r(rc)
                #         else:
                #             ghost_negative -= -coeff * r(rc)
                #     continue
                for (f11, f12), v1 in cem_tensor1.items():
                    for (f21, f22), v2 in cem_tensor2.items():
                        first_result = next(iter((g(f11) * g(f21)).to_rc_graph_ring_element())).resize(n)
                        result += v1 * v2 * r(first_result.squash_product(f12.resize(n)).squash_product(f22.resize(n)))
                # rc_prd = (g.from_dict(cem_dict1) * g.from_dict(cem_dict2)).to_rc_graph_ring_element()
                # # assert len(rc_prd) == 1, f"Multiple terms in product for {perm1}, {perm2} at {rc1.perm}, {rc2.perm}: {rc_prd}"
                # # assert next(iter(rc_prd.values())) == 1, f"Coefficient not 1 in product for {perm1}, {perm2} at {rc1.perm}, {rc2.perm}: {rc_prd}"
                # result += rc_prd
        #assert ghost_negative.almosteq(r.zero), f"Negative ghost terms remaining for {perm1}, {perm2}: {ghost_negative}, {ghost_positive}"
        #assert ghost_positive.almosteq(r.zero), f"Positive ghost terms remaining for {perm1}, {perm2}: {ghost_positive}, {ghost_negative}"
        #result += ghost_positive - ghost_negative
        for rc, coeff in result.items():
            assert coeff == prd.get(rc.perm, 0), f"Failure for {perm1}, {perm2} at {rc.perm} with coeff {coeff} vs {prd.get(rc.perm, 0)}"
        print("Success for", perm1, perm2)

        
        
