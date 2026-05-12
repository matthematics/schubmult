from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.qy_rc_graph_ring import _canonical_rc
from schubmult.rings.combinatorial.forest_rc_ring import _canonical_rc as forcanonical_rc
from schubmult.utils.perm_utils import weak_compositions, artin_sequences

def verify_key_polynomial_weight_formulas(n):
    print("Verifying key polynomial weight formulas...")
    comps = artin_sequences(n - 1)
    the_canonical_words = {}
    KeyDual = FreeAlgebra(KeyBasis)
    r = RCGraphRing()
    for comp in comps:
        the_monomial = r.monomial(*comp)
        for rc in the_monomial.keys():
            if rc.extremal_weight not in the_canonical_words:
                the_canonical_words[rc.extremal_weight] = rc.to_highest_weight()[0].perm_word
            elif rc.to_highest_weight()[0].perm_word < the_canonical_words[rc.extremal_weight]:
                the_canonical_words[rc.extremal_weight] = rc.to_highest_weight()[0].perm_word
    for comp in comps:
        the_monomial = r.monomial(*comp)
        result = 0
        for rc, coeff in the_monomial.items():
            if rc.to_highest_weight()[0].perm_word == the_canonical_words[rc.extremal_weight]:
                result += coeff * KeyDual(*rc.extremal_weight)
        assert result.change_basis(WordBasis).almosteq(FA(*comp)), f"Failed for {comp} with result {result.change_basis(WordBasis)}"
    print("Key polynomial weight formulas verified successfully!")

def verify_slide_polynomial_weight_formulas(n):
    print("Verifying slide polynomial weight formulas...")
    comps = artin_sequences(n - 1)
    the_canonical_words = {}
    FSlideDual = FreeAlgebra(FundamentalSlideBasis)
    r = RCGraphRing()
    for comp in comps:
        the_monomial = r.monomial(*comp)
        for rc in the_monomial.keys():
            can_rc = _canonical_rc(rc)
            if can_rc.length_vector not in the_canonical_words:
                the_canonical_words[can_rc.length_vector] = can_rc.perm_word
            elif can_rc.perm_word < the_canonical_words[can_rc.length_vector]:
                the_canonical_words[can_rc.length_vector] = can_rc.perm_word
    for comp in comps:
        the_monomial = r.monomial(*comp)
        result = 0
        for rc, coeff in the_monomial.items():
            can_rc = _canonical_rc(rc)
            if rc.perm_word == the_canonical_words[can_rc.length_vector]:
                result += coeff * FSlideDual(*can_rc.length_vector)
        assert result.change_basis(WordBasis).almosteq(FA(*comp)), f"Failed for {comp} with result {result.change_basis(WordBasis)}"
    print("Slide polynomial weight formulas verified successfully!")

def verify_forest_polynomial_weight_formulas(n):
    print("Verifying forest polynomial weight formulas...")
    comps = artin_sequences(n - 1)
    the_canonical_words = {}
    ForestDual = FreeAlgebra(ForestBasis)
    r = RCGraphRing()
    for comp in comps:
        the_monomial = r.monomial(*comp)
        for rc in the_monomial.keys():
            can_rc = forcanonical_rc(rc)
            if can_rc.forest_weight not in the_canonical_words:
                the_canonical_words[can_rc.forest_weight] = can_rc.perm_word
            elif can_rc.perm_word < the_canonical_words[can_rc.forest_weight]:
                the_canonical_words[can_rc.forest_weight] = can_rc.perm_word
    for comp in comps:
        the_monomial = r.monomial(*comp)
        result = 0
        for rc, coeff in the_monomial.items():
            can_rc = forcanonical_rc(rc)
            if can_rc.perm_word == the_canonical_words[can_rc.forest_weight]:
                result += coeff * ForestDual(*can_rc.forest_weight)
        assert result.change_basis(WordBasis).almosteq(FA(*comp)), f"Failed for {comp} with result {result.change_basis(WordBasis)}"
    print("Forest polynomial weight formulas verified successfully!")

def quasi_rcs(perm, length, weight=None):
    return {rc for rc in RCGraph.all_rc_graphs(perm, length,weight=weight) if rc.is_quasi_yamanouchi}

def verify_slide_schub_prod_coeff_formula(n):
    import itertools
    print("Verifying slide polynomial/Schubert polynomial coefficient formulas...")
    comps = artin_sequences(n - 1)
    perms = Permutation.all_permutations(n)
    FSlideDual = FreeAlgebra(FundamentalSlideBasis)
    r = RCGraphRing()
    for w in perms:
        for rc in {_canonical_rc(rc0) for rc0 in RCGraph.all_rc_graphs(w)}:
            c = rc.length_vector
            cprd = FSlideDual(*c).coproduct()
            cprd_schub = ASx(w, len(rc)).coproduct()
            for ((u, _), (v, _)), schub_coeff in cprd_schub.items():
                the_sum = 0
                the_prod = Sx(u) * Sx(v)
                for pram, scoeff in the_prod.items():
                    the_sum += scoeff * len(quasi_rcs(pram, len(rc), c))
                coeff_magic = 0
                for (a, b), coeff in cprd.items():
                    coeff_magic += coeff * len(quasi_rcs(u, len(rc), a))*len(quasi_rcs(v, len(rc), b))
                assert the_sum == coeff_magic    
    print("Fancy verified successfully!")

if __name__ == "__main__":
    import sys
    n = int(sys.argv[1])
    #max_degree = int(sys.argv[2])

    # verify_slide_polynomial_weight_formulas(n)
    # verify_forest_polynomial_weight_formulas(n)
    # verify_key_polynomial_weight_formulas(n)
    verify_slide_schub_prod_coeff_formula(n)
    