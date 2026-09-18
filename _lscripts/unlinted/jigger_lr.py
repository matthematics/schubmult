from schubmult import *
import itertools

def all_squashes(rc1, n):
    perms = Permutation.all_permutations(n)
    tensors = set()
    for perm in perms:
        for rc3 in RCGraph.all_rc_graphs(perm, len(rc1), weight=rc1.length_vector):
            cems = RCGraph.full_CEM(rc3.perm, n, partition=tuple(range(n-1, 0, -1)))
            if rc1 in cems:
                tensors.update(cems[rc1].keys())
    return tensors

def coprod(rc):
    g = GrassTensorAlgebra()
    r = RCGraphRing()
    if isinstance(rc, RCGraph):
        if rc.perm.inv == 0:
            return g(()) @ g(())
        perm = rc.perm
        n = len(rc)
        N = len(rc.perm.trimcode)
        cprd_try = (g@g).zero
        perms = Permutation.all_permutations(len(rc.perm))
        for perm1, perm2 in itertools.product(perms, repeat=2):
            if perm1.inv + perm2.inv != perm.inv:
                continue
            if len(perm1.trimcode) > N or len(perm2.trimcode) > N:
                continue

            for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n), RCGraph.all_rc_graphs(perm2, n)):
                if rc.length_vector != tuple([ a + b for a, b in zip(rc1.length_vector, rc2.length_vector)]):
                    continue
                squash1 = all_squashes(rc1, n)
                squash2 = all_squashes(rc2, n)
                for sq1, sq2 in itertools.product(squash1, squash2):
                    elem  = (g(sq1) * g(sq2)).to_rc_graph_ring_element()
                    if elem.resize(n).almosteq(r(rc)):
                        cprd_try += g(sq1) @ g(sq2)
        return cprd_try
    result = (g@g).zero
    for rc1, coeff1 in rc.items():
        result += coeff1 * coprod(rc1)
    return result

if __name__ == "__main__":
    import sys
    

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    g = GrassTensorAlgebra()
    r = RCGraphRing()
    for perm1, perm2 in itertools.product(perms, repeat=2):
        if perm1.inv == 0 or perm2.inv == 0:
            continue
        for rc1, rc2 in itertools.product(RCGraph.all_rc_graphs(perm1, n), RCGraph.all_rc_graphs(perm2, n)):
            cp1 = coprod(rc1)
            cp2 = coprod(rc2)
            cp12 = coprod(r(rc1)*r(rc2))
            assert cp12.almosteq(cp1 * cp2), f"Failure for {rc1}, {rc2} with coprod {cp12}, expected {cp1 * cp2}\n{cp12 - cp1*cp2}"
    # for perm in perms:
    #     for rc in RCGraph.all_rc_graphs(perm, n):
    #         base_coprod = coprod(r(rc))
    #         left_coprod = (r@r@r).zero
    #         right_coprod = (r@r@r).zero
    #         for (rc1, rc2), coeff in base_coprod.items():
    #             left_coprod += coeff * coprod(r(rc1)) @ r(rc2)
    #             right_coprod += coeff * r(rc1) @ coprod(r(rc2))
    #         assert left_coprod.almosteq(right_coprod), f"Failure for {rc} with left coprod {left_coprod}, right coprod {right_coprod}\nDifference: {left_coprod - right_coprod}"
    #     print("Success for", perm)
        # rc = RCGraph.principal_rc(perm, n)
        
        # bumble = (ASx@ASx).zero
        # for (rc1,rc2), coeff in cprd_try.items():
        #     bumble += r(rc1).resize(n).to_free_algebra_element() @ r(rc2).resize(n).to_free_algebra_element()
        # real_cprd = ASx(perm, n).coproduct()
        # assert bumble.almosteq(real_cprd), f"Failure for {perm} with cprd {bumble}, expected {real_cprd}\n{bumble-real_cprd}"
        print("Success for", perm1, perm2)