from schubmult import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import *
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
    """Grove polynomial of ``comp`` built by inverting the omega insertion.

    We enumerate the compatible set-valued labelings ``kappa`` of the indexed
    forest ``F = weak_composition_to_indfor(comp)`` (the grove definition), and
    realize each labeling as a WCGraph by writing down its (word, compatible
    sequence) pair *explicitly* -- the published inverse of the set-valued omega
    insertion -- and feeding it to ``WCGraph.from_word_compatible``.

    Concretely:

    * The canonical left-binary-search labeling ``P`` of ``F`` and the decreasing
      labeling ``Q`` of the principal reduced RC graph form the omega pair for
      ``F``. The map ``Gamma = omega_reduced_word_from_labelings`` reads the
      reduced word ``W`` off ``(P, Q)`` -- this is the inverse of the insertion,
      so ``W`` is determined by ``F`` (not chosen at random). Node ``v`` sits at
      word position ``len(W) - Q(v)`` and carries the reduced letter
      ``ell(v) = W[len(W) - Q(v)]``.
    * A labeling ``kappa`` places the letter ``ell(v)`` into every row
      ``r in kappa(v)``. Reading the resulting ``(row, letter)`` pairs in
      ``(row, -letter)`` order gives a weakly-increasing compatible sequence and
      a word whose rows are strictly decreasing, i.e. exactly the data
      ``WCGraph.from_word_compatible`` consumes. The multiset of compatible
      values is the disjoint union of the ``kappa(v)``, so the WCGraph monomial
      equals ``x^kappa``.

    Each WCGraph contributes ``beta**(|kappa| - |F|)`` times its monomial; beta is
    the degree ``-1`` homogenizer recording extra labels beyond one per node.
    """
    from itertools import combinations

    from schubmult.combinatorics.rc_graph import RCGraph

    forest = weak_composition_to_indfor(comp)

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

    # Omega pair (P, Q) for F, read off the principal reduced RC graph.
    principal = RCGraph.principal_rc(uncode(comp), length)
    p_labeling, q_labeling = principal.omega_invariant

    # Gamma (the published inverse of the insertion) reconstructs the reduced word
    # from (P, Q); reversing matches the forward convention of perm_word.
    word = list(reversed(omega_reduced_word_from_labelings(p_labeling, q_labeling)))

    # The reduced letter carried by node v is W at v's Q-position.
    letter_of = {
        node.index: word[len(word) - q_labeling(node.index)]
        for node in p_labeling.forest.inorder_traversal
    }

    result = 0
    wc_set = set()
    for labeling in all_labelings:
        # Explicit (word, compatible sequence): node v's letter into each row of kappa(v).
        entries = [(row, letter_of[index]) for index, label_set in labeling.items() for row in label_set]
        entries.sort(key=lambda pair: (pair[0], -pair[1]))
        compat_seq = [row for row, _ in entries]
        wc_word = [letter for _, letter in entries]
        wc = WCGraph.from_word_compatible(wc_word, compat_seq, length=length)
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
        
        #grove_try = grove_rc_try(comp, Gx._beta, n - 1)
        # forest_dict = grove_as_forest_dict(comp, beta=Gx._beta, length=n - 1)
        # forest_result = Forest.from_dict(forest_dict)
        # forest_rcs = {k for k in RCGraph.all_rc_graphs(uncode(comp), n - 1) if k.is_forest_rc}
        # fatpants_sum = 0
        # for dinkbat in forest_rcs:
        #     new_comp = dinkbat.forest_weight
        #     grove_result = grove_polynomial(new_comp, Sx.genset, Gx._beta)
        #     grove_try = sum([wc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for wc in WCGraph.grove_wcs(new_comp, n - 1, dinkbat)])
        #     diff =  (grove_result - grove_try).expand() 
        #     assert diff == 0, f"Mismatch for {new_comp}: {grove_result} != {grove_try}\n{dinkbat=}\nDiff: {diff}"
        #     print(f"Puncho! {dinkbat.is_principal=}")
        #     fatpants_sum += grove_try
        # print("Hofer gonk")

        wc_set = WCGraph.all_wc_graphs(uncode(comp), n - 1)
        # for dinkbat in forest_rcs:
        #     wc_set -= WCGraph.grove_wcs(dinkbat.forest_weight, n - 1, dinkbat)
        fatpants_sum = 0
        # while wc_set:
        new_forest_rcs = {wc for wc in wc_set if wc.forest_weight == wc.length_vector}
        for dinkbat in new_forest_rcs:
            gang = WCGraph.grove_wcs(dinkbat.forest_weight, n - 1, dinkbat)
            #assert gang.issubset(wc_set), f"Mismatch for {dinkbat.forest_weight}: {gang.difference(wc_set)}"
            #wc_set -= gang
            grover = sum([wc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for wc in gang])
            assert (grover - (Gx._beta ** (len(dinkbat.perm_word) - sum(comp))) * grove_polynomial(dinkbat.forest_weight, Sx.genset, Gx._beta)).expand() == 0, f"Mismatch for {dinkbat.forest_weight}: {grover} != {(Gx._beta ** (len(gang.pop().perm_word) - sum(dinkbat.forest_weight))) * grove_polynomial(dinkbat.forest_weight, Sx.genset, Gx._beta)}"
            print("Profie")
            fatpants_sum += grover
            #sum([wc.polyvalue(Sx.genset, beta=Gx._beta, prop_beta=True) for wc in gang])
        assert (fatpants_sum - grothendieck_poly(uncode(comp), Sx.genset, ZeroGeneratingSet(), Gx._beta)).expand() == 0, f"Mismatch for {comp}: {fatpants_sum} != {grothendieck_poly(uncode(comp), Sx.genset, ZeroGeneratingSet(), Gx._beta)}\nDiff: {(fatpants_sum - grothendieck_poly(uncode(comp), Sx.genset, ZeroGeneratingSet(), Gx._beta)).expand()}"
        print("Potato piston")
        # grothy = grothendieck_poly(uncode(comp), Sx.genset, ZeroGeneratingSet(), Gx._beta)
        # assert (fatpants_sum - grothy).expand() == 0, f"Mismatch for {comp}: {fatpants_sum} != {grothy}\nDiff: {(fatpants_sum - grothy).expand()}"
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