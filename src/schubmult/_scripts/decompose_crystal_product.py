import itertools

from schubmult import CrystalGraph, CrystalGraphTensor
from schubmult import RCGraph
from schubmult import RCGraphRing

from schubmult import *

if __name__ == "__main__":
    import sys

    from schubmult.utils.perm_utils import elem_sym_perms
    from sympy import pretty_print

    n = int(sys.argv[1])
    rc_ring = RCGraphRing()
    perms = Permutation.all_permutations(n)
    # for perm1, perm2 in itertools.product(perms, perms):
    #     graphs1 = RCGraph.all_rc_graphs(perm1, n-1)
    #     graphs2 = RCGraph.all_rc_graphs(perm2, n-1)
    #     dct = {}
    #     hw_check = set()
    #     for rc1, rc2 in itertools.product(graphs1, graphs2):
    #         cprd = CrystalGraphTensor(rc1, rc2).to_highest_weight()[0]
    #         if cprd in hw_check:
    #             continue
    #         hw_check.add(cprd)
    #         prod1 = rc_ring(rc1)*rc_ring(rc2)
    #         tring = rc_ring @ rc_ring
    #         raise_elem = tring.zero
    #         print(f"Coproduct of {perm1}, {perm2} highest weights")
    #         for tens in cprd.full_crystal:
    #             raise_elem = rc_ring.zero
    #             rc01, rc02 = tens.factors[0], tens.factors[1]
    #             prod2 = rc_ring(rc01)*rc_ring(rc02)
    #             for rc, coeff in prod2.items():
    #                 hw = CrystalGraphTensor(rc01, rc02).to_highest_weight()[0]
    #                 raise_elem += coeff * rc_ring(hw.factors[0]) * rc_ring(hw.factors[1])
    #             df = raise_elem-prod1
    #             print(f"{df=}")
    #             try:
    #                 assert all(v == 0 for v in df.values()), f"Mismatch in coproduct for {perm1}, {perm2} at highest weight {cprd} {raise_elem=} {prod1=}"

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


    hw_elems = {}
    lr = {}
    tring = rc_ring @ rc_ring
    w0_rc = RCGraph.principal_rc(Permutation.w0(n), n-1)
    for perm in perms:

        #for rc in RCGraph.all_rc_graphs(perm, n-1):
        hw = set()
        hw2 = set()
        results = {}
        product = Sx(perm)*Sx(w0_rc.perm)
        for rc_p in RCGraph.all_rc_graphs(perm,n-1):
            prod = rc_ring(rc_p) * rc_ring(w0_rc)
            thw = CrystalGraphTensor(w0_rc,rc_p).to_highest_weight()[0]
            hw2.add(thw.factors[1])
            for rc, coeff in prod.items():
                # pretty_print(rc)
                # print(f"{len(rc)=} {coeff=}")
                tp = CrystalGraphTensor(RCGraphCut(rc, n-1, len(rc)), RCGraphCut(rc, 0, n-1))
                rc_hw = tp.to_highest_weight()[0]
                frc = RCGraph([*rc_hw.factors[0].under_rc[:n-1],*rc_hw.factors[1].under_rc[n-1:]])
                hw.add(rc_hw.factors[1].base_graph)

                # results[full_rc_hw] = results.get(full_rc_hw, set())
                # results[full_rc_hw].add(rc_p)
        print(f"Results for perm {perm}:")
        print("hw")
        for rc in hw:
            pretty_print(rc)
        print("hw2")
        for rc in hw2:
            pretty_print(rc)

        # for rc_cry, rc_result in results.items():
        #     #print(f"{rc_result.inverse_crystal.is_lowest_weight=}")
        #     pretty_print(rc_cry)
        #     pretty_print(rc_cry.crystal_weight)
            # for rc_p in rc_result:
            #     pretty_print(rc_p)
            #     wt = tuple(a+b for a,b in zip(rc_p.length_vector, w0_rc.length_vector))
            #     print(wt)
            #     print(uncode(wt))
            #     print(f"{product.get(uncode(wt), 0)=}")
        print(f"{len(product)=}")
        print(f"{len(results)=}")
        print(F"{product=}")
        input()
        # try:
        #     assert all(v == 0 for v in (cprd - our_cprd).values()), f"Mismatch in coproduct for perm {perm}"
        # except AssertionError as e:
        #     print(e)
        #     print(cprd - our_cprd)
        #     input()
        # for shift in range(len(prin_rc)):
        #     rc = RCGraph([*([()]* shift), *[tuple([a + len(prin_rc) for a in row]) for row in prin_rc], *([()]*(len(prin_rc) - shift))]).transpose().resize(2*len(prin_rc))
        #     i = len(prin_rc)
        #     hw = CrystalGraphTensor(RCGraphCut(rc, 0, i), RCGraphCut(rc, i, len(rc))).to_highest_weight()
        #     hw_rc = RCGraph([*hw[0].factors[0].under_rc.rowrange(0, i), *hw[0].factors[1].under_rc[i:]])
        #     rc1 = hw[0].factors[0].under_rc.vertical_cut(i)[0]

        #     rc2 = hw[0].factors[1].under_rc.rowrange(i)
        #     if (rc1, rc2) in the_set:
        #         continue
        #     lr[perm] = lr.get(perm, tring.zero) + tring((rc1, rc2))
        #     the_set.add((rc1, rc2))
        #     hw_elems[(rc1, rc2)] = hw_elems.get((rc1, rc2), set())
        #     hw_elems[(rc1, rc2)].add(hw_rc)

        #     hw2 = CrystalGraphTensor(*rc.vertical_cut(i)).to_highest_weight()
        #     assert hw2[0].factors == (rc1, rc2)
    # for (rc1, rc2), elems in hw_elems.items():
    #     print(f"Highest weights for")
    #     pretty_print(CrystalGraphTensor(rc1, rc2))
    #     print(lr[(~rc1.perm, ~rc2.perm)])
    for perm in lr:
        elem = lr[perm]
        print(f"Perm {perm=}")
        pretty_print(lr[perm])
        print(ASx(~perm, n-1).coproduct())
            # for elem in elems:
            #     assert elem.is_valid

