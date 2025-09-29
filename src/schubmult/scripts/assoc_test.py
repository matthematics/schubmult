# Traceback (most recent call last):
#   File "/home/matthematics/schubmult/src/schubmult/scripts/assoc_test.py", line 163, in <module>
#     diff = hom(g1 * g2) - hom(g1) * hom(g2)
#                ~~~^~~~
#   File "/home/matthematics/schubmult/src/schubmult/rings/rc_graph_module.py", line 656, in __mul__
#     return self.prod_with_rc(other)
#            ^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/home/matthematics/schubmult/src/schubmult/rings/rc_graph_module.py", line 1631, in prod_with_rc
#     new_buildup_module += RCGraphModule(dict.fromkeys(rc.right_zero_act(), coeff))
#                                                       ^^^^^^^^^^^^^^^^^^^
#   File "/home/matthematics/schubmult/src/schubmult/rings/rc_graph_module.py", line 1129, in right_zero_act
#     assert len(rc_set) == len(up_perms), f"{rc_set=}, {len(up_perms)=}, {self=} {up_perms=}"
#            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# AssertionError: rc_set={RCGraph((1,), (2,), ()), RCGraph((1,), (3,), ())}, len(up_perms)=3, self=RCGraph((1,), (2,)) up_perms=ASx(((1, 3, 4, 2), 3)) + ASx(((2, 1, 4, 3), 3)) + ASx(((2, 3, 1), 3))
# we need mulling 0 to give us 2 3

def hom3(rc):
    from schubmult import FA, ASx, SchubertBasis
    from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
    if isinstance(rc, RCGraph):
        #return (ASx@FA)(r,rc.length_vector()))
        ring = ASx @ ASx
        #first = FA(*rc.length_vector()).change_basis(SchubertBasis).get()
        return ring.ext_multiply(ASx(rc.perm,len(rc)),ASx(rc.perm,len(rc)))
    ret = 0
    for rc0, coeff in rc.items():
        ret += coeff * hom(rc0)
    return ret

def hom_rc(elem):
    from schubmult import FA, ASx, SchubertBasis
    from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
    ret = 0
    for key, val in elem.items():
        ret+=val * RCGraphModule(dict.fromkeys({RCGraph.principal_rc(*key)},1))
    return ret


def hom(rc):
    from schubmult import FA, ASx, SchubertBasis
    from schubmult.rings import TensorRing
    from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
    ring = TensorRing(ASx, RCGraphModule())
    
    if isinstance(rc, RCGraph):
        #return (ASx@FA)(r,rc.length_vector()))
        # print(f"{rc.length_vector()} {tuple(rc)=}")    
        #first = FA(*rc.length_vector()).change_basis(SchubertBasis).get()
        return ASx(rc.perm,len(rc))
        #return ret
    ret = 0
    for rc0, coeff in rc.items():
        # print(f"{rc=}")
        # print(f"{w0=}")
        # print(f"{rc0=}")
        ret += coeff * hom(rc0)
    return ret

# def hom(rc):
#     from schubmult import FA, ASx, SchubertBasis
#     from schubmult.rings import TensorRing
#     from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
#     ring = TensorRing(ASx, RCGraphModule())
    
#     if isinstance(rc, RCGraph):
#         #return (ASx@FA)(r,rc.length_vector()))
#         # print(f"{rc.length_vector()} {tuple(rc)=}")    
#         #first = FA(*rc.length_vector()).change_basis(SchubertBasis).get()
#         return ASx(rc.perm,len(rc))
#         #return ret
#     ret = 0
#     for rc0, coeff in rc.items():
#         # print(f"{rc=}")
#         # print(f"{w0=}")
#         # print(f"{rc0=}")
#         ret += coeff * hom(rc0)
#     return ret


def fa_hom(rc):
    from schubmult import FA, ASx, SchubertBasis
    from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
    if isinstance(rc, RCGraph):
        #return (ASx@FA)(r,rc.length_vector()))
        ring = FA @ (ASx@ASx)
        first = FA(*rc.length_vector()).change_basis(SchubertBasis)
        second = FA(*rc.length_vector()).change_basis(SchubertBasis).coproduct()
        return ring.ext_multiply(first, second)
    ret = 0
    for rc0, coeff in rc.items():
        ret += coeff * hom_cop(rc0)
    return ret


def hom_cop(rc):
    from schubmult import FA, ASx, SchubertBasis
    from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule
    if isinstance(rc, RCGraph):
        #return (ASx@FA)(r,rc.length_vector()))
        ring = ASx @ (ASx@ASx)
        first = FA(*rc.length_vector()).change_basis(SchubertBasis)
        second = FA(*rc.length_vector()).change_basis(SchubertBasis).coproduct()
        return ring.ext_multiply(first, second)
    ret = 0
    for rc0, coeff in rc.items():
        ret += coeff * hom_cop(rc0)
    return ret



if __name__ == "__main__":
    # test module functionality
    
    import sys

    from schubmult import FA, ASx, Permutation, SchubertBasis, WordBasis, uncode
    from schubmult.abc import x
    from schubmult.rings import MonomialBasis, PolynomialAlgebra, SchubertPolyBasis
    from schubmult.rings.rc_graph_module import RCGraph, RCGraphModule, all_fa_degree
    from schubmult.utils.perm_utils import artin_sequences

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    perms = Permutation.all_permutations(n)
    modfull = 0
    len1 = 4
    dct = {}
    deg = 6

    # if we coprod a Schub this will act right

    

    # for perm in perms:
    #     AG = PolynomialAlgebra(SchubertPolyBasis(len(perm.trimcode)))
    #     bob = AG(perm).coproduct()
    #     mod = 0
    #     print(perm)
    #     for (perm1, perm2), v in bob.items():
    #         rc1 = RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(perm1[0], perm1[1]),1))
    #         rc2 = RCGraphModule(dict.fromkeys(RCGraph.all_rc_graphs(perm2[0], perm2[1]),1))
    #         print(f"{perm1=} {perm2=} {v=}")
    #         print(v * (rc1 * rc2))
    #         mod = v * (rc1 * rc2)
    #         for k, v in mod.items():
    #             assert v == 0 or k.perm == perm, f"{k=}, {perm=}, {v=}, {k.perm=}"
    #     #print(mod)
        
    #for deg in range((n*(n-1))//2+1):
    # for seq in all_fa_degree(deg, len1):
    #     g1 = FA(*seq[:3]) * RCGraph()
    #     g2 = FA(*seq[3:]) * RCGraph()
    #     print(g1 * g2)
    #     print(FA(*seq).change_basis(SchubertBasis))

    # for perm in perms:
    #     print(f"{perm=}")
    #     mod = 0
    #     wrd = ASx(perm, len1).change_basis(WordBasis)
    #     for w, v in wrd.items():
    #         if w in dct:
    #             mod += hom(v * dct[w])
    #     print(mod)
    for len1 in range(1,n-1):
        for perm in perms:
            if len(perm.trimcode) > len1:
                continue
            graphs1 = RCGraph.all_rc_graphs(perm, len1)
            for len2 in range(1,n-1):
                for perm2 in perms:
                    if len(perm2.trimcode) > len2:
                        continue
                    graphs2 = RCGraph.all_rc_graphs(perm2, len2)
                    # g1 = graphs1
                    # g2 = graphs2
                    
                    for g1 in graphs1:
                        for g2 in graphs2:
                            diff = hom(g1 * g2) - hom(g1) * hom(g2)
                            print(f"Testing {perm} in {len1}, {perm2} in {len2}")
                            print("g1")
                            print(hom(g1))
                            print("g2")
                            print(hom(g2))
                            print("g1*g2")
                            print(g1*g2)
                            print(hom(g1)*hom(g2))
                            print("the_hom")
                            print(hom(g1*g2))
                            # print(f"{(g1*g2).value_dict=}")
                            # print(f"{dict(hom(g1)*hom(g2))=}")
                            assert all(v == 0 for k, v in diff.items()), f"{tuple(diff.items())=}"
                            print("Success")
                    # elem = ASx(perm1, len1)
                    # elem2 = ASx(perm2, len2)
                    # print(f"{hom3(elem).value_dict} {hom3(elem2).value_dict}")
                    # print(f"{(hom3(elem)*hom3(elem2)).value_dict} - {hom3(elem*elem2).value_dict}")
                    # diff = hom3(elem)*hom3(elem2) - hom3(elem*elem2)
                    # print(diff)
                    # assert all(v == 0 for v in diff.value_dict.values())
                    # print("Success")