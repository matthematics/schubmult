from schubmult import *
from schubmult.combinatorics.indexed_forests import *
from schubmult.rings.combinatorial.forest_rc_ring import ForestRCGraphRing
from schubmult.rings.polynomial_algebra import *
from schubmult.symbolic.common_polys import *
from schubmult.symbolic.poly.variables import *
from schubmult.rings.polynomial_algebra import *
from schubmult.utils._mul_utils import add_perm_dict
from schubmult.utils.tuple_utils import pad_tuple

import itertools

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
    """GrovePoly polynomial of ``comp`` built by inverting the omega insertion.

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

def _grove_it_up(comp, bw, n):
    grove_asgroth = GrovePoly(*comp).change_basis(GrothendieckPolyBasis)
    # groth = bw.full_groth_elem(uncode(comp), n, 1)
    grove = 0
    for (perm, length), coeff in grove_asgroth.items():
        grove += coeff * bw.full_groth_elem(perm, n, 1)
        # rc_key = bw.key_to_wc_graph(key).resize(len(comp))
        # if rc_key.grove_weight == tuple(comp) or rc_key.perm != uncode(comp):
        # grove += coeff * bw(key)
    #grove = bw.from_dict({k: v for k, v in grove.items() if (bw.key_to_wc_graph(k).resize(len(comp)).grove_weight == tuple(comp) or bw.key_to_wc_graph(k).perm != uncode(comp))})
    return grove

if __name__ == "__main__":
    import sys

    bw = BoundedWCFactorAlgebra()
    n = int(sys.argv[1])
    perms = Permutation.all_permutations(n)
    comps = [tuple(perm.pad_code(n -  1)) for perm in perms]
    the_poles = {}
    for comp1, comp2 in itertools.product(comps, repeat=2):

        if comp1 not in the_poles:
            grove1 = _grove_it_up(comp1, bw, n)
            the_poles[comp1] = grove1
        else:
            grove1 = the_poles[comp1]
        if comp2 not in the_poles:
            grove2 = _grove_it_up(comp2, bw, n)
            the_poles[comp2] = grove2
        else:
            grove2 = the_poles[comp2]

        producto = (grove1 * grove2).to_wc_graph_ring_element().resize(n - 1)

        real_prod = GrovePoly(*comp1) * GrovePoly(*comp2)

        checko_prod = 0
        for wc, v in producto.items():
            if wc.forest_weight == wc.length_vector:
                checko_prod += v * GrovePoly(*wc.forest_weight)

        assert real_prod.almosteq(checko_prod), f"Failed for {comp1} * {comp2}: {real_prod-checko_prod=}"
        print("Pantoopa fatcough")