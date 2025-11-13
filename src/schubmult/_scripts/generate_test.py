from schubmult import RCGraphRing
from sympy import pretty_print

# IS THIS ASSOCIATIVE?
# need associativity
rc_ring = RCGraphRing()

def hom3(rc):
    from schubmult import FA, ASx, SchubertBasis
    #from schubmult import RCGraph, rc_ring.from_dict
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
    #from schubmult import RCGraph, rc_ring.from_dict
    ring = RCGraphRing()
    ret = 0
    for key, val in elem.items():
        ret+=val * rc_ring.from_dict(dict.fromkeys({RCGraph.principal_rc(*key)},1))
    return ret


def hom(rc):
    from schubmult import TensorRing

    from schubmult import FA, ASx, SchubertBasis
    #from schubmult import RCGraph, rc_ring.from_dict
    ring = TensorRing(ASx@FA, rc_ring)

    if isinstance(rc, RCGraph):
        #return (ASx@FA)(r,rc.length_vector()))
        # print(f"{rc.length_vector()} {tuple(rc)=}")
        #first = FA(*rc.length_vector()).change_basis(SchubertBasis).get()
        return (ASx@FA)(((rc.perm,len(rc)),rc.length_vector()))
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
#     from schubmult import TensorRing
#     from schubmult import RCGraph, rc_ring.from_dict
#     ring = TensorRing(ASx, rc_ring.from_dict())

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
    # from schubmult import RCGraph, rc_ring.from_dict
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
    # from schubmult import RCGraph, rc_ring.from_dict
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


def single_rc(a):
    if a == 0:
        return RCGraph(((),))
    return RCGraph([tuple(range(a, 0, -1))])

def RC(*seq):
    from schubmult import RCGraph
    res = rc_ring.from_dict({RCGraph(): 1})
    for a in seq:
        res = res * single_rc(a)
    return res


def csym_rc(*weight):
    buildup = []
    for i, a in enumerate(weight):
        if a == 0:
            buildup.append(())
        else:
            buildup.append(tuple(range(a+i, i, -1)))
    return rc_ring.from_dict({RCGraph(buildup): 1})

if __name__ == "__main__":
    # test module functionality

    import sys

    from schubmult import x
    from schubmult import MonomialBasis, PolynomialAlgebra, SchubertPolyBasis
    from schubmult import RCGraph
    from schubmult.utils.perm_utils import artin_sequences

    from schubmult import FA, ASx, Permutation, SchubertBasis, WordBasis, uncode

    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    perms = Permutation.all_permutations(n)
    modfull = 0
    len1 = 4
    dct = {}
    deg = 6
    rc_ring = RCGraphRing()

    # if we coprod a Schub this will act right
    # print(csym_rc(1)*csym_rc(2))
    # print(csym_rc(1)*csym_rc(2) - csym_rc(1,2) )
    # exit()

    # for perm in perms:
    #     AG = PolynomialAlgebra(SchubertPolyBasis(len(perm.trimcode)))
    #     bob = AG(perm).coproduct()
    #     mod = 0
    #     print(perm)
    #     for (perm1, perm2), v in bob.items():
    #         rc1 = rc_ring.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm1[0], perm1[1]),1))
    #         rc2 = rc_ring.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm2[0], perm2[1]),1))
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
    # rc1 = RCGraph([(1,)])
    # rc2 = RCGraph([(3,),(2,)])
    # print(rc1)
    # print(rc2)
    # print(rc1 * rc2)
    # print(rc2 * rc1)
    # exit()
    # rc_ring = RCGraphRing()
    # for perm in perms:
    #     elem = ASx(perm).change_basis(WordBasis)
    #     mod = 0
    #     rc = next(iter(RCGraph.all_rc_graphs(perm)))
    #     print(rc)
    #     print(elem)
    #     for w, v in elem.items():
    #         elem2 = rc_ring.from_dict({RCGraph(): v})
    #         index = 0
    #         for a in w:
    #             elem2 =  elem2*csym_rc(*list((rc.length_vector()[index:index+a])))
    #             index += a
    #         mod += elem2
    #     #mod = RC(*perm.trimcode)
    #     print(f"{perm.trimcode}")
    #     print(mod)
    # exit()
    from schubmult import Sx

    def elem_func(p, k, *args):
        if p == k and p == 0:
            return 1
        if k < p or p < 0:
            return 0
        return rc_ring(RCGraph.one_row(p))

    longest = uncode(list(range(n-1,0,-1)))
    for perm in perms:
        pretty_print(Sx(perm).cem_rep(mumu=longest, elem_func=elem_func))
