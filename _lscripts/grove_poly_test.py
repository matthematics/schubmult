from schubmult import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.poly.schub_poly import *
from schubmult.symbolic.poly.variables import *
from schubmult.utils._mul_utils import add_perm_dict
from schubmult.utils.tuple_utils import pad_tuple

f = ForestRCGraphRing()

T = ThompsonAlgebra()

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

# def grove_as_thompson(comp, beta, length):
#     forest = weak_composition_to_indfor(comp)
#     desc = forest.trim_descents[0]
#     is_left_child = forest.is_left_child(desc)
#     trimmed_forest = forest.trim_descent(desc)
#     return T((desc,)) * ( grove_as_thompson(trimmed_forest.code, beta, length) 
def act_on_grove_dict(dct, beta, length):
    result = {}
    for comp, coeff in dct.items():
        forest = weak_composition_to_indfor(comp)
        if len(forest.trim_descents) == 0:
            result[comp] = result.get(comp, 0) + coeff
            continue
        desc = forest.trim_descents[0]
        trimmed_forest = forest.trim_descent(desc)
        trimmed_comp = trimmed_forest.code
        result[trimmed_comp] = result.get(trimmed_comp, 0) + coeff
        result[comp] = result.get(comp, 0) + coeff * beta * T((-desc,))
    return result

def grove_as_forest_dict(comp, beta, length):
    thmp = grove_as_thompson(comp, beta, length)
    print(f"Thompson result for {comp}: {thmp}")    
    result = {}
    for comp, coeff in thmp.items():
        comp0 = tuple([c for c in comp if c > 0])
        forest_comp = pad_tuple(indexed_forest_from_trimming_word(comp0).code, length)
        result[forest_comp] = result.get(forest_comp, 0) + coeff
    return result


def groth_as_rc(perm, beta, length):
    dct = WCGraph.groth_to_schub(perm, beta=beta)
    result = f.zero
    for perm, coeff in dct.items():
        # for rc in RCGraph.all_rc_graphs(perm, length):
        #     if rc.is_forest_rc:
        result += f.from_dict(dict.fromkeys(RCGraph.all_rc_graphs(perm, length), coeff), snap=True)
    return result

def _forest_weight(wc):
    if wc.perm.inv == len(wc.perm_word):
        return RCGraph(wc).forest_weight
    return wc._snap_reduced().forest_weight

def grove_rc_try(comp, beta, length):
    """Grove polynomial of ``comp`` built by going from the indexed forest to
    WCGraphs.

    We enumerate the compatible set-valued labelings of the indexed forest
    ``F = weak_composition_to_indfor(comp)`` (Definition of grove polynomials),
    and realize each labeling as a WCGraph on the principal reduced word of
    ``comp``. The forest node -> word-position map comes from the ``omega``
    ``Q``-labeling of that reduced RC graph, and each labeling's node sets become
    the set-sequence of the WCGraph. Each WCGraph contributes
    ``beta**(|kappa| - |F|)`` times its monomial, i.e. beta is the degree ``-1``
    homogenizer recording extra labels beyond one per node.
    """
    from itertools import combinations

    from schubmult.combinatorics.rc_graph import RCGraph

    forest = weak_composition_to_indfor(comp)
    size = len(forest)

    def labelings_below(node, min_value):
        results = []
        for subset_size in range(1, node.rho - min_value + 2):
            for subset in combinations(range(min_value, node.rho + 1), subset_size):
                top = subset[-1]
                left_opts = labelings_below(node.left, top) if node.left is not None else [{}]
                right_opts = labelings_below(node.right, top + 1) if node.right is not None else [{}]
                for left in left_opts:
                    for right in right_opts:
                        labeling = {node.index: tuple(subset)}
                        labeling.update(left)
                        labeling.update(right)
                        results.append(labeling)
        return results

    all_labelings = [{}]
    for root in forest._roots:
        root_labelings = labelings_below(root, 1)
        all_labelings = [{**base, **choice} for base in all_labelings for choice in root_labelings]

    # Principal reduced RC graph for comp gives the reduced word and, via its
    # omega Q-labeling, the forest node -> word position correspondence.
    principal = next(
        rc
        for rc in RCGraph.all_rc_graphs(uncode(comp), length)
        if rc.forest_weight == tuple(comp) and uncode(comp) == rc.perm
    )
    word = list(principal.perm_word)
    _, q_labeling = principal.omega_invariant
    position_of = {
        node.index: len(word) - q_labeling(node.index)
        for node in principal.forest_invariant.forest.inorder_traversal
    }

    result = 0
    wc_set = set()
    for labeling in all_labelings:
        set_sequence = [None] * len(word)
        total_labels = 0
        for index, label_set in labeling.items():
            set_sequence[position_of[index]] = set(label_set)
        wc = WCGraph.from_reduced_compatible_set_sequence(word, set_sequence, length=length)
        assert wc.is_valid
        wc_set.add(wc)
    for wc in wc_set:
        result += (beta ** (len(wc.perm_word) - sum(comp))) * wc.polyvalue(Sx.genset)
    return result

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [perm.pad_code(n -  1) for perm in perms]
    for comp in comps:
        grove_result = grove_polynomial(comp, Sx.genset, Gx._beta)
        grove_try = grove_rc_try(comp, Gx._beta, n - 1)
        # forest_dict = grove_as_forest_dict(comp, beta=Gx._beta, length=n - 1)
        # forest_result = Forest.from_dict(forest_dict)
        diff =  (grove_result - grove_try).expand() 
        assert diff == 0, f"Mismatch for {comp}: {grove_result} != {grove_try}\nDiff: {diff}"
        print("Hofer gonk")
    # for w in perms:
    #     if w.inv == 0:
    #         continue
        # groth_result = grothendieck_poly(w, Sx.genset, ZeroGeneratingSet(), Gx._beta)
        # spinach = groth_as_rc(w, beta=Gx._beta, length=n - 1)
        # assert (spinach.polyvalue(Sx.genset) - groth_result).expand() == 0, f"Mismatch for {w}: {spinach.polyvalue(Sx.genset)} != {groth_result}"
        # print(spinach)
        # print("This is certainly a baked potato")
        # comps = set([rc.forest_weight for rc in spinach])
        # result = 0
        # for comp in comps:
        #     word = []
        #     forest = weak_composition_to_indfor(comp)
        #     while forest.trim_descents:
        #         desc = forest.trim_descents[0]
        #         word = [(desc, forest.is_left_child(desc)), *word]
        #         forest = forest.trim_descent(desc)
        #     result += _g_grove_extractor(spinach, indexes=word, beta=Gx._beta) * grove_polynomial(comp, Sx.genset, Gx._beta)
        
        # assert (result - groth_result).expand() == 0, f"Mismatch for {w}: {result=}, {groth_result=}\ndiff={(result - groth_result).expand()}"
        #print(f"Success for {w}")

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