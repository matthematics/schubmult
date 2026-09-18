from schubmult import *
from schubmult.rings.polynomial_algebra import *

def cauchy_dual_key(n):
    """Demazure atoms?"""
    perms = Permutation.all_permutations(n)
    # seen = {}
    build_polys = {}
    build_second_polys = {}
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            # if rc.extremal_weight != RCGraph.principal_rc(perm, n - 1).length_vector:
            #     continue
            complementary_weight = [n - 1 - i - rc.length_vector[i] for i in range(n - 1)]
            dict_key = (tuple(rc.perm.pad_code(n - 1)),rc.to_highest_weight()[0])
            build_polys[dict_key] = build_polys.get(dict_key, 0) + PA(*complementary_weight)
            # dict_right_key = rc.to_highest_weight()[0]
            # build_second_polys[dict_right_key] = build_second_polys.get(dict_right_key, 0) + PA(*rc.length_vector)
    for wtt in sorted(build_polys.keys(), key=lambda f: f[1].extremal_weight):
        poly = build_polys[wtt]
        
        print(f"{wtt[0]},{wtt[1].extremal_weight}: {poly.change_basis(KeyPolyBasis)}") 
    # for keykey, poly1 in build_second_polys.items():
    #     print(f"foundbat    {keykey.extremal_weight}: {poly1.change_basis(KeyPolyBasis)}")

def cauchy_dual_forest(n):
    """???"""
    from schubmult.rings.combinatorial.forest_rc_ring import _canonical_rc
    perms = Permutation.all_permutations(n)
    # seen = {}
    build_polys = {}
    are_these_the_same = {}
    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n - 1):
            # if rc.extremal_weight != RCGraph.principal_rc(perm, n - 1).length_vector:
            #     continue
            complementary_weight = [n - 1 - i - rc.length_vector[i] for i in range(n - 1)]
            dict_key = (tuple(rc.perm.pad_code(n - 1)),_canonical_rc(rc))
            build_polys[dict_key] = build_polys.get(dict_key, 0) + PA(*complementary_weight)
            #stump_key
    for wtt in sorted(build_polys.keys(), key=lambda f: f[1].forest_weight):
        poly = build_polys[wtt]
        print(f"{wtt[0]},{wtt[1].forest_weight}: {poly}") 


if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    cauchy_dual_key(n)
    # cauchy_dual_forest(n)