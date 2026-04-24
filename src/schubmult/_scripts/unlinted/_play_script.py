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

def convert_tensor(tensor):
        ret = (ASx@FA@ASx@FA).zero
        for (rc1, rc2), coeff in tensor.items():
            ret += coeff * ASx(rc1.perm, len(rc1)) @ FA(*rc1.length_vector) @ ASx(rc2.perm, len(rc2)) @ FA(*rc2.length_vector)
        return ret

def schub_weight_equal(tensor1, tensor2):
    
    
    return convert_tensor(tensor1).almosteq(convert_tensor(tensor2))

neg_pats = [Permutation([1, 4, 3, 2]), Permutation([4, 1, 3, 2]), Permutation([3, 1, 4, 2])]
r = RCGraphRing()
def rc_to_ss(rc1):
    ret = EE.zero
    if isinstance(rc1, RCGraph):
        schuby = r(rc1).to_free_algebra_element().change_basis(ElementaryBasis)
        wordy = FA(*rc1.length_vector).change_basis(ElementaryBasis)
        
        for elem_comp, coeff in schuby.items():
            if coeff == 0:
                continue
            if wordy.get(elem_comp, 0) != 0:
                ret += coeff * EE(*elem_comp)
    else:
        for rc, coeff in rc1.items():
            ret += coeff * rc_to_ss(rc)
    return ret

def full_coprod(perm, weight, length):
    cem = RCGraph.full_CEM(perm, length)
    result = 0
    for rc, cem_dict in cem.items():
        if rc.length_vector != weight:
            continue
        for rc_tup, coeff in cem_dict.items():
            if coeff == 0:
                continue
            if tuple([a + b for a, b in zip(rc_tup[0].length_vector, rc_tup[1].length_vector)]) == weight:
                yield (rc_tup, coeff)

ring = GrassTensorAlgebra()
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
            grass_elem += ring.from_dict(cem_dict)
        elif rc.perm != perm:
            grass_elem += ring.from_dict(cem_dict)
        else:
            rc_base, rc_grass = rc.resize(n).squash_decomp()
            if rc_grass.perm.inv != 0:
                grass_elem += ring.from_dict({(*k, rc_grass): v for k, v in RCGraph.full_CEM(rc_base.perm, n - 1)[rc_base.resize(n - 1)].items()})
            else:
                grass_elem += ring.from_dict(cem_dict)
    return grass_elem

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


def elem_key_to_dual_algebra(key):
    index = 0
    build_tuple = []
    length = key.size
    for rc22 in key:
        if len(rc22) == revpar[index]:
            build_tuple.append(rc2.perm.inv)
            index += 1
        elif index >= partition[0]:
            break
        else:
            # diff = len(rc2) - 1 - revpar[index]
            # build_tuple.extend(0 for j in range(diff))
            # index += diff
            # build_tuple.append(rc2.perm.inv)
            while revpar[index] < len(rc2):
                build_tuple.append(0)                            
                index += 1
            build_tuple.append(rc2.perm.inv)
            index += 1
    if index < partition[0]:
        build_tuple.extend([0] * len(revpar[index:]))
    the_tup = tuple(pad_tuple(build_tuple, partition[0]))
    the_tup = the_tup[:length - 1] + tuple(sorted(the_tup[length - 1:]))
    dual_elem2 = EE(the_tup, n - 1)
    return dual_elem2

if __name__ == "__main__":
    from schubmult.combinatorics.crystal_graph import CrystalGraphTensor
    from schubmult.abc import *
    from schubmult.rings.polynomial_algebra import *
    from sympy import pretty_print
    import sys
    n = int(sys.argv[1])
    extra = 0
    perms = Permutation.all_permutations(n)

    r = BoundedRCFactorAlgebra()

    for perm in perms:
        spit = r.schub_elem(perm, len(perm.trimcode), partition=tuple((~(perm.mul_dominant())).trimcode))
        if any(v < 0 for v in spit.values()):
            assert any(perm.has_pattern(neg_pat) for neg_pat in [[3,1,4,2],[4,1,3,2],[1,4,3,2]]), f"Failure for {perm}, got {spit}"
        # else:
        #     assert not any(perm.has_pattern(neg_pat) for neg_pat in [[3,1,4,2],[4,1,3,2],[1,4,3,2]]), f"Failure for {perm}, got {spit}"
    exit(0)
    # def _is_elem_sym(rc):
    #     return rc.perm.descents() == {len(rc) - 1} and set(rc.trimcode).issubset({0, 1})
    # EE = FreeAlgebra(ElementaryBasis)
    # for p2 in range(2, n):
    #     for k in range(p2, n):
    #         p1 = k - 1
    #         for rc1_key, rc2_key in itertools.product(r.elem_sym(p1, k - 1).keys(), r.elem_sym(p2, k).keys()):
    #             rc1, rc2 = r.key_to_rc_graph(rc1_key), r.key_to_rc_graph(rc2_key)
    #             squash = rc1.resize(len(rc2)).squash_product(rc2)
    #             if len(squash.perm.trimcode) <= k and squash.is_highest_weight:
    #                 rc_base, rc_grass = squash.resize(k).squash_decomp()
    #                 assert rc_grass.inv == rc2.inv, f"Failure for {p1}, {p2}, got {rc1}, {rc2}, squash {squash}, base {rc_base}, grass {rc_grass}"
    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     length = len(perm.trimcode)
    #     #schub_elem = r.schub_elem(perm, length)
    #     partition = tuple((~(perm.strict_mul_dominant())).trimcode)
    #     if partition[-1] > 1:
    #         partition += (tuple(range(partition[-1] - 1, 0, -1)))

    #     cem = RCGraph.full_CEM(perm, length, partition=partition)
    #     revpar = tuple(reversed(partition))
    #     print(revpar)
    #     length2 = len(revpar)
    #     the_dual_elem = ASx(perm, length).change_basis(ElementaryBasis)
    #     dct = {}
    #     for rc, cem_dict in cem.items():
            
    #         for key, coeff in cem_dict.items():
    #             build_tuple = []
    #             index = 0
    #             for rc22 in key:
    #                 if len(rc22) == revpar[index]:
    #                     build_tuple.append(rc2.perm.inv)
    #                     index += 1
    #                 elif index >= partition[0]:
    #                     break
    #                 else:
    #                     # diff = len(rc2) - 1 - revpar[index]
    #                     # build_tuple.extend(0 for j in range(diff))
    #                     # index += diff
    #                     # build_tuple.append(rc2.perm.inv)
    #                     while revpar[index] < len(rc2):
    #                         build_tuple.append(0)                            
    #                         index += 1
    #                     build_tuple.append(rc2.perm.inv)
    #                     index += 1
    #             if index < partition[0]:
    #                 build_tuple.extend([0] * len(revpar[index:]))
    #             the_tup = tuple(pad_tuple(build_tuple, partition[0]))
    #             the_tup = the_tup[:length - 1] + tuple(sorted(the_tup[length - 1:]))
    #             check_perm = (rc.perm* (~(perm.mul_dominant())))
                
                
                
    #             #assert (the_dual_elem.get((check_perm, len(the_tup)), 0) != 0 and len(cem_dict) == 1) or (the_dual_elem.get((check_perm, len(the_tup)), 0) == 0 and len(cem_dict) > 1), f"Failure for {perm}, got {key}, which corresponds to {rc.perm=} {rc=} with tuple {the_tup}, {coeff=}, {the_dual_elem=} {check_perm=}"
    #             assert (the_dual_elem.get((the_tup, length), 0) > 0 and rc.perm == perm) or (rc.perm != perm), f"Failure for {perm}, got {key}, which corresponds to {rc.perm=} {rc=} with tuple {the_tup}, {coeff=}, {the_dual_elem=} {check_perm=}"
    #             #if rc.perm == perm:
    #         # dct[(the_tup, length)] = dct.get((the_tup, length), 0) + coeff
    #         # the_dual_elem2 = EE.from_dict(dct)
    #     #assert the_dual_elem.almosteq(the_dual_elem2), f"Failure for {perm}, got {key}, which corresponds to {rc.perm=} {rc=} with tuple {the_tup}, {coeff=}, {the_dual_elem=} {the_dual_elem2=}"
    #         #((~other_perm) * perm).is_dominant, f"Failure for {perm}, got {key}, which corresponds to {other_perm} with tuple {the_tup}, {coeff=}, {(~other_perm) * perm=}"
    #         #assert tuple([rc.perm.inv for rc in key]) == pad_tuple(tuple(perm.trimcode), length), f"Failure for {perm}, got {key}"
    
    

    # # grass_tensor_elems = {}
    # # hw_grass = {}
    r = RCGraphRing()
    g = BoundedRCFactorAlgebra()
    def elem_sym_perm(p, k):
        return uncode(([0] * (k - p)) + [1]*p)

    # for perm in perms:
    #     if perm.inv == 0:
    #         continue
    #     #max_desc = max(perm.descents()) + 1
    #     elem = g.schub_elem(perm, len(perm.trimcode) + 1)
    #     # grass blast
    #     # new_elem = g.zero
    #     # for key, coeff in elem.items():
    #     #     grass_key = key[-1]
    #     #     new_elem += coeff * g(g.make_key((*key[:-1],RCGraph(), key.size)) * g.schub_elem(grass_key.perm, key.size)
    #     # assert new_elem.almosteq(elem), f"Failure for {perm}, got {elem}, got {new_elem} from grass blast"
    #     pretty_print(elem)

    for perm, p2 in itertools.product(perms, range(1, n)):
        # if p1 >= p2:
        #     continue
        for k2 in range(p2, n):
            # if k1 >= k2:
            #     continue
            #perm1 = elem_sym_perm(p1, k1)
            perm2 = elem_sym_perm(p2, k2)
            length = max(k2, len(perm.trimcode))
            elem1 = g.schub_elem(perm, length)
            elem2 = g.schub_elem(perm2, length)
            prod = elem1 * elem2
            prod2 = elem2 * elem1
            assert prod2.almosteq(prod), f"Failure for {perm1} and {perm2}, got {prod}, got {prod2} from reversed product"
            
            # for key, coeff in prod.items():
            #     if coeff == 0:
            #         continue
            #     check_perm = uncode(([0] * (key.size - p1 - p2)) + [1]*(p1 + p2))
            #     assert check_perm in prod, f"Failure for {perm1} and {perm2}, got {prod}, missing {check_perm}"