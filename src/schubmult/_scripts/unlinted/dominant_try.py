import itertools

from schubmult.ck_ring import CoxeterKnuthRing
from schubmult import CrystalGraph, CrystalGraphTensor


class RCGraphCut(CrystalGraph):

    def __init__(self, under_rc, start, end):
        self.under_rc = under_rc
        self.start = start
        self.end = end

    @property
    def crystal_weight(self):
        return self.base_graph.crystal_weight

    def crystal_length(self):
        return self.end - self.start

    def epsilon(self, i):
        return self.base_graph.epsilon(i)

    def phi(self, i):
        return self.base_graph.phi(i)

    def raising_operator(self, index):
        index += self.start
        if index < 1 or index >= self.end:
            return None
        new_val = self.under_rc.raising_operator(index)
        if new_val is None:
            return None
        return RCGraphCut(new_val, self.start, self.end)

    def lowering_operator(self, index):
        index += self.start
        if index < 1 or index >= self.end:
            return None
        new_val = self.under_rc.lowering_operator(index)
        if new_val is None:
            return None
        return RCGraphCut(new_val, self.start, self.end)

    @property
    def base_graph(self):
        bg = self.under_rc.rowrange(self.start)
        if self.end - self.start < len(bg):
            bg = bg.vertical_cut(self.end - self.start)[0]
        return bg

if __name__ == "__main__":
    import sys

    from schubmult import RCGraph
    from schubmult import RCGraphRing
    from schubmult.utils.perm_utils import elem_sym_perms
    from sympy import pretty_print

    from schubmult import *

    n = int(sys.argv[1])
    rc_ring = RCGraphRing()
    perms = Permutation.all_permutations(n)
    # cd = []
    # for i in range(2 * (n - 1), 0, -2):
    #     cd += [i]
    # TRANSPORT CRYSTAL
    dominant_graphs = {RCGraph.principal_rc(perm.minimal_dominant_above(), n-1) for perm in perms}
    #perms2n = {perm for perm in Permutation.all_permutations(2 * n - 1) if perm.bruhat_leq(uncode(cd))}
    for dom in dominant_graphs:
        if dom.perm.inv == 0:
            continue
        for perm in perms:
            if perm.inv == 0:
                continue
            product = Sx(perm) * Sx(dom.perm)
            schubs = ASx(dom.perm, n - 1) * ASx(perm, n - 1)
            zob = Sx.zero
            print(f"{zob=}")
            for (perm1, _), coeff in schubs.items():
                zob += coeff * Sx(perm1*(~uncode([*((0,)*(n-1)), *perm.trimcode])))
            results = Sx.zero
            nhw = 0
            for perm2, coeff in zob.items():
                dct = dict(Sx(perm2).coproduct(*list(range(1,n-1))))
                print(f"{dct=}")
                print(f"{dom.perm,perm=}")
                nhw += coeff * dct.get((dom.perm,Permutation([])), 0)
            for rc in RCGraph.all_rc_graphs(perm, n - 1):

                rcs = dom.prod_with_rc(rc)
                for graph in rcs:
                    zop = 0
                    tp = CrystalGraphTensor(RCGraphCut(graph, 0, n-1), RCGraphCut(rc, n-1, len(rc)))
                    for tp_hw in tp.all_highest_weights():
                        if tp_hw.factors[0].base_graph == dom:
                            zop += 1
                    assert nhw == zop, f"Mismatch in dominant transport test for {perm}, {dom.perm} {nhw=} {len(tp.all_highest_weights())=}"
                    print("Success")
