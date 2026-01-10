from schubmult import *

def fa_elem(rc):
    FA = FreeAlgebra(WordBasis)
    ASx = FreeAlgebra(SchubertBasis)
    if isinstance(rc, RCGraph):        
        return (ASx@FA)(((rc.perm, len(rc)), rc.length_vector))
    ring = ASx @ FA
    return ring.from_dict({next(iter(fa_elem(rc).keys())): coeff for rc, coeff in rc.items()})

def killifbad(fa):
    spitball = fa.ring.zero
    bang = False
    for ((perm,ln), word), coeff in fa.items():
        if len(RCGraph.all_rc_graphs(perm, ln, weight=word)) == 0:
            bang = True
            continue
        spitball += coeff * fa.ring(((perm,ln),word))
    if bang:
        print(f"{spitball=} {fa=}")
    return spitball

def rc_elem(fa_elem):
    ring = RCGraphRing()
    sm = ring.zero
    for ((perm, ln), word), coeff in fa_elem.items():
        sm += ring.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, ln, weight=word), coeff))
    return sm

def da_schub_(perm, length):
    ring = RCGraphRing()
    return ring.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, length), 1))

def da_schub(fa_elem, length):
    if isinstance(fa_elem, Permutation):
        return da_schub_(fa_elem, length)
    ring = RCGraphRing()
    sm = ring.zero
    for (perm, ln), coeff in fa_elem.items():
        sm += coeff * da_schub_(perm, ln)
    return sm

if __name__ == "__main__":
    import sys
    import itertools
    from symengine import S
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    rc_set = set()
    for perm in perms:
        rc_set.update(RCGraph.all_rc_graphs(perm, n - 1))

    fa_indexes = set()
    for rc in rc_set:
        fa_indexes.update(fa_elem(rc).keys())


    # ring = RCGraphRing()
    for rc1, rc2 in itertools.product(rc_set, repeat=2):
        print(f"Testing rc1: {rc1.perm}, rc2: {rc2.perm}")
        fa1 = fa_elem(rc1)
        fa2 = fa_elem(rc2)
        print(f"In da stink {rc1*rc2=}")
        print(f"In da fish {fa1 * fa2=}")

    # print("GOOD")
    # tring = TensorRing(ASx, FA)
    # for perm1, perm2 in itertools.product(perms, repeat=2):
    #     da1 = da_schub(perm1, n - 1)
    #     da2 = da_schub(perm2, n - 1)
    #     prod_da = da1 * da2
    #     conv_prod_da = sum([coeff * da_schub(*p) for p, coeff in (ASx(perm1, n-1) * ASx(perm2, n-1)).items()])
    #     assert all(v == S.Zero for v in (prod_da - conv_prod_da).values()), f"Mismatch:\nfa: {prod_da}\ntr: {conv_prod_da}"
    print("GOOD afsopin")