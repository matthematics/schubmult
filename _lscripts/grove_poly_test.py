from schubmult import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.symbolic.poly.schub_poly import *
from schubmult.symbolic.poly.variables import *


f = ForestRCGraphRing()

def _g_operator(rc, index, beta):
    """The operator G_i on rc graphs, which is the composition of the quasi-shift and trim-descent operations."""
    return f(rc).forest_trim(index) - beta * f(rc).quasi_shift(index)

def _g_grove_extractor(rc_elem, indexes, beta):
    """The operator G_i on rc graphs, which is the composition of the quasi-shift and trim-descent operations."""
    result = f.zero
    if len(indexes) > 0:
        desc, is_left_child = indexes[-1]      
        for rc, coeff in rc_elem.items():
            forest = weak_composition_to_indfor(rc.forest_weight)
            if desc not in forest.trim_descents:
                continue
            if is_left_child:
                result += coeff * (_g_operator(rc, desc, beta) + beta * f(rc).quasi_shift(desc))
            else:
                result += coeff * (_g_operator(rc, desc, beta) + beta * f(rc).quasi_shift(desc + 1))
        if len(indexes) > 1:
            return _g_grove_extractor(result, indexes=indexes[:-1], beta=beta)
    else:
        result = rc_elem
    result_values = [v for k, v in result.items() if k.perm.inv == 0]
    return sum(result_values)

def groth_as_rc(perm, beta, length):
    dct = WCGraph.groth_to_schub(perm, beta=beta)
    result = f.zero
    for perm, coeff in dct.items():
        # for rc in RCGraph.all_rc_graphs(perm, length):
        #     if rc.is_forest_rc:
        result += f.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, length), coeff), snap=True)
    return result

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    
    for w in perms:
        if w.inv == 0:
            continue
        groth_result = grothendieck_poly(w, Sx.genset, ZeroGeneratingSet(), Gx._beta)
        spinach = groth_as_rc(w, beta=Gx._beta, length=n - 1)
        assert (spinach.polyvalue(Sx.genset) - groth_result).expand() == 0, f"Mismatch for {w}: {spinach.polyvalue(Sx.genset)} != {groth_result}"
        print(spinach)
        print("This is certainly a baked potato")
        comps = set([rc.forest_weight for rc in spinach])
        result = 0
        for comp in comps:
            word = []
            forest = weak_composition_to_indfor(comp)
            while forest.trim_descents:
                desc = forest.trim_descents[0]
                word = [(desc, forest.is_left_child(desc)), *word]
                forest = forest.trim_descent(desc)
            result += _g_grove_extractor(spinach, indexes=word, beta=Gx._beta) * grove_polynomial(comp, Sx.genset, Gx._beta)
        
        assert (result - groth_result).expand() == 0, f"Mismatch for {w}: {result=}, {groth_result=}\ndiff={(result - groth_result).expand()}"
        print(f"Success for {w}")

        # fordict = {}
        # for wc in WCGraph.all_wc_graphs(w, n - 1):
        #     word = []
        #     compat = []
        #     perm_last = Permutation([])
        #     letter_last = 0
        #     append_val = 1
        #     for a in wc.perm_word:
        #         new_perm = perm_last @ Permutation.ref_product(a)
        #         if new_perm != perm_last:
        #             word.append(a)
        #             perm_last = new_perm
        #             if a >= letter_last:
        #                 append_val += 1
        #             compat.append(append_val)
        #             letter_last = a
        #     rc = RCGraph.from_reduced_compatible(word, compat)
        #     fordict[rc.forest_invariant] = fordict.get(rc.forest_invariant, 0) + wc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True)
        # testo = {}
        # for invar, val in fordict.items():
        #     cd = invar.forest.code
        #     if cd in testo:
        #         assert (testo[cd] - val).expand() == 0, f"Mismatch for invariant {invar} code {cd}: {testo[cd]} != {val}"
        #     else:
        #         testo[cd] = val
        # print("happy")