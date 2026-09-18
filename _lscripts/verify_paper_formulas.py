from schubmult import *
from schubmult.rings.free_algebra import *
from schubmult.rings.combinatorial.qy_rc_graph_ring import _canonical_rc
from schubmult.rings.combinatorial.forest_rc_ring import _canonical_rc as forcanonical_rc
from schubmult.utils.perm_utils import weak_compositions, artin_sequences

def verify_key_polynomial_lr_rule(N, M, r, R):
    print("Verifying key polynomial Littlewood-Richardson rule...")
    for n in range(r, N + 1):
        comps1 = weak_compositions(n, R)
        for m in range(r, M + 1):
            comps2 = weak_compositions(m, R)
            KeyDual = FreeAlgebra(KeyBasis)
            r = RCGraphRing()
            for comp1 in comps1:
                key1 = KeyDual(*comp1)
                for comp2 in comps2:
                    key2 = KeyDual(*comp2)
                    prd = key1 * key2
                    for extwt, coeff in prd.items():
                        #result += coeff * KeyDual(*rc.extremal_weight)
                        C = next(iter([ac for ac in r.monomial(*extwt).keys() if ac.extremal_weight == extwt]))
                        result = 0
                        for rc in C.full_crystal:
                            top, bottom = rc.vertical_cut(len(comp1))
                            
                            if top.extremal_weight == comp1 and bottom.extremal_weight == comp2 and top.is_highest_weight and bottom.is_highest_weight:
                                result += 1
                        assert result == coeff, f"Failed for {comp1} and {comp2} with result {result} and coeff {coeff}"
                        print(f"Hot bang potato {comp1} {comp2} {extwt} with result {result} and coeff {coeff}")
            #assert result.change_basis(WordBasis).almosteq(FA(*comp1)*FA(*comp2)), f"Failed for {comp1} and {comp2} with result {result.change_basis(WordBasis)}"
    print("Key polynomial Littlewood-Richardson rule verified successfully!")

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

def verify_forest_polynomial_lr_rule(N, M):
    print("Verifying forest polynomial Littlewood-Richardson rule...")
    for n in range(1, N + 1):
        comps1 = artin_sequences(n - 1)
        for m in range(1, M + 1):
            comps2 = artin_sequences(m - 1)
            ForestDual = FreeAlgebra(ForestBasis)
            r = RCGraphRing()
            for comp1 in comps1:
                forest1 = ForestDual(*comp1)
                for comp2 in comps2:
                    forest2 = ForestDual(*comp2)
                    prd = forest1 * forest2
                    for extwt, coeff in prd.items():
                        #result += coeff * ForestDual(*rc.forest_weight)
                        monom = r.monomial(*extwt)
                        C = next(iter([ac for ac in monom.keys() if ac.forest_weight == extwt]))
                        result = 0
                        for rc in {rcc for rcc in RCGraph.all_rc_graphs(C.perm, len(C)) if rcc.forest_invariant == C.forest_invariant}:
                            if rc.forest_weight != extwt:
                                continue
                            top, bottom = rc.vertical_cut(len(comp1))
                            
                            if top.forest_weight == comp1 and bottom.forest_weight == comp2 and top.length_vector == comp1 and bottom.length_vector == comp2:
                                result += 1
                        assert result == coeff, f"Failed for {comp1} and {comp2} with result {result} and coeff {coeff}"
                        print(f"Hot bang potato {comp1} {comp2} with result {result} and coeff {coeff}")
            #assert result.change_basis(WordBasis).almosteq(FA(*comp1)*FA(*comp2)), f"Failed for {comp1} and {comp2} with result {result.change_basis(WordBasis)}"
    print("Forest polynomial Littlewood-Richardson rule verified successfully!")

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
    m = int(sys.argv[2])
    r = int(sys.argv[3])
    R = int(sys.argv[4])
    #max_degree = int(sys.argv[2])

    # verify_slide_polynomial_weight_formulas(n)
    # verify_forest_polynomial_weight_formulas(n)
    # verify_key_polynomial_weight_formulas(n)
    # verify_slide_schub_prod_coeff_formula(n)
    verify_key_polynomial_lr_rule(n, m, r, R)
    # verify_forest_polynomial_lr_rule(n, m)
    