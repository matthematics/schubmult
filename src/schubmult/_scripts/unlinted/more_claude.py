from schubmult import *
from sympy import pretty_print
"""
Refined understanding of divided differences and crystal structure.

Key correction: divdiff_desc(i) ≠ 0 requires THREE conditions:
1. i is a descent of rc.perm
2. The simple root α_i appears in row i of the RC-graph
3. Deleting α_i yields the i-highest weight

The result is the i-highest weight (root deleted), and divdiff(i) 
includes the full i-string starting from this highest weight.
"""

def clarify_divdiff_structure():
    """
    Clarify the structure of divided differences on RC-graphs.
    
    CORRECT understanding:
    
    1. divdiff_desc(rc, i) returns non-zero IFF:
       a) i ∈ Des(rc.perm)  (i is a descent)
       b) α_i is in row i of rc  (the simple root is present)
       c) Deleting α_i yields the i-highest weight
    
    2. When non-zero, divdiff_desc(rc, i) returns:
       - The RC-graph with α_i deleted from row i
       - This is the i-highest weight of the i-string
    
    3. divdiff(rc, i) returns:
       - The full i-string starting from divdiff_desc(rc, i)
       - This is: hw + f_i(hw) + f_i²(hw) + ... 
       - Where hw = divdiff_desc(rc, i)
       - And f_i is the crystal lowering operator
    
    4. The i-string condition:
       - divdiff_desc gives the HIGHEST weight of the string
       - The string continues via lowering operators
       - If deletion doesn't give highest weight, divdiff_desc = 0
    """
    print("\n" + "="*80)
    print("CORRECT STRUCTURE OF DIVIDED DIFFERENCES")
    print("="*80)
    
    print("\ndivdiff_desc(rc, i) ≠ 0  IFF:")
    print("  1. i ∈ Des(rc.perm)")
    print("  2. α_i appears in row i of rc")
    print("  3. Deleting α_i yields the i-highest weight")
    print()
    print("When these hold:")
    print("  divdiff_desc(rc, i) = rc with α_i deleted from row i")
    print()
    print("divdiff(rc, i) includes:")
    print("  - hw = divdiff_desc(rc, i)")
    print("  - f_i(hw), f_i²(hw), ..., f_i^k(hw)")
    print("  - Until f_i returns None")
    print()
    print("Key insight: divdiff_desc tests the i-string structure!")
    print("For the SAME permutation, different RC-graphs can have different divdiff!")
    print()


def when_is_divdiff_zero(n):
    """
    Analyze when divdiff(rc, i) = 0 for ALL RC-graphs, not just principal ones.
    
    This happens when:
    1. rc has no descent at i, OR
    2. rc has descent at i but α_i is not in row i, OR
    3. Deleting α_i doesn't yield i-highest weight
    
    CORRECTED: divdiff_desc returns a SET, not a single element.
    Empty set means divdiff_desc = 0.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"WHEN IS divdiff(rc, i) = 0? (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    zero_cases = {
        'no_descent': [],
        'no_root': [],
    }
    
    # Generate ALL RC-graphs using the library function
    all_rcs = []
    for perm in perms:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    print(f"\nAnalyzing {len(all_rcs)} total RC-graphs")
    print()
    
    for rc in all_rcs:
        for i in range(1, n):
            diff = rc_ring(rc).divdiff(i)
            
            if len(diff) == 0:
                # Classify why it's zero
                if i not in rc.perm.descents():
                    zero_cases['no_descent'].append((rc, i))
                else:
                    # i is a descent, check divdiff_desc
                    # divdiff_desc returns a SET - empty means it's 0
                    desc_result = rc.divdiff_desc(i)
                    
                    if len(desc_result) == 0:  # Empty set
                        # α_i not in row i, or doesn't yield hw
                        zero_cases['no_root'].append((rc, i))
                    else:
                        # divdiff_desc is non-empty but full divdiff is empty
                        # This is unexpected - means the i-string has no elements?
                        print(f"  ⚠ INTERESTING: RC-graph with perm {rc.perm.trimcode}, i={i}")
                        print(f"    divdiff = 0 (empty) but divdiff_desc ≠ 0 (has {len(desc_result)} elements)")
                        print(f"    This might mean the i-string starting from divdiff_desc is empty?")
                        print(f"    divdiff_desc result:")
                        for desc_rc in desc_result:
                            print(f"      perm {desc_rc.perm.trimcode}, weight {desc_rc.crystal_weight}")
                        pretty_print(rc)
                        print()
    
    print(f"\nCases where divdiff(rc, i) = 0:")
    print()
    print(f"1. No descent at i: {len(zero_cases['no_descent'])} cases")
    print(f"2. Descent but divdiff_desc also 0: {len(zero_cases['no_root'])} cases")
    print()
    
    # Show ALL examples
    print("\n" + "="*60)
    print("ALL 'no descent' cases:")
    print("="*60)
    for rc, i in zero_cases['no_descent']:
        print(f"\nRC {rc.perm.trimcode}: i={i}")
        print(f"  Descents: {rc.perm.descents()}")
        print(f"  i={i} is NOT a descent (so divdiff must be 0)")
        print(f"  Weight: {rc.crystal_weight}")
        pretty_print(rc)
        print()
    
    print("\n" + "="*60)
    print("ALL 'descent but divdiff_desc = 0' cases:")
    print("="*60)
    for rc, i in zero_cases['no_root']:
        print(f"\nRC {rc.perm.trimcode}: i={i}")
        print(f"  Descents: {rc.perm.descents()}")
        print(f"  i={i} IS a descent, but:")
        print(f"    divdiff_desc({i}) = empty set")
        print(f"    Meaning: α_i not in row i, or deletion doesn't give i-hw")
        print(f"  Weight: {rc.crystal_weight}")
        pretty_print(rc)
        print()
    
    return zero_cases


def find_useful_zero_divdiff_pairs(n):
    """
    Find pairs (rc, i) where divdiff(rc, i) = 0 AND this is useful
    for extending products via Leibniz rule.
    
    KEY: Check ALL RC-graphs, not just principal ones!
    
    CORRECTED: divdiff_desc returns a SET.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"USEFUL ZERO DIVDIFF PAIRS FOR EXTENSION (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant]
    rc_ring = RCGraphRing()
    
    # Generate ALL RC-graphs using the library function
    all_rcs = []
    for perm in perms:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    print(f"\nAnalyzing {len(all_rcs)} total RC-graphs")
    print(f"Looking for pairs (rc, i) where:")
    print(f"  - divdiff(rc, i) = 0")
    print(f"  - ∃ dominant dom with divdiff(dom, i) ≠ 0")
    print()
    
    useful_pairs = []
    
    for rc in all_rcs:
        for i in range(1, n):
            diff_rc = rc_ring(rc).divdiff(i)
            
            if len(diff_rc) == 0:
                # Check if this helps extend products
                # Useful if there exists a dominant dom such that ∂_i(dom) ≠ 0
                
                can_extend = False
                for dom_perm in dom_perms:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    
                    if len(diff_dom) > 0:
                        can_extend = True
                        break
                
                if can_extend:
                    useful_pairs.append((rc, i))
    
    print(f"\nFound {len(useful_pairs)} useful pairs")
    print()
    
    # Group by permutation
    by_perm = {}
    for rc, i in useful_pairs:
        if rc.perm not in by_perm:
            by_perm[rc.perm] = {}
        if i not in by_perm[rc.perm]:
            by_perm[rc.perm][i] = []
        by_perm[rc.perm][i].append(rc)
    
    print("="*60)
    print("COMPLETE LISTING OF ALL USEFUL PAIRS")
    print("="*60)
    print()
    
    for perm in sorted(by_perm.keys(), key=lambda p: (p.inv, p.trimcode)):
        print(f"\nPermutation {perm.trimcode}:")
        print(f"  Inversion count: {perm.inv}")
        print(f"  Descents: {perm.descents()}")
        print(f"  Is dominant: {perm.is_dominant}")
        
        principal = RCGraph.principal_rc(perm, n-1)
        
        for i in sorted(by_perm[perm].keys()):
            rcs = by_perm[perm][i]
            
            print(f"\n  Index i={i}:")
            print(f"    {len(rcs)} RC-graph(s) with divdiff(rc, {i}) = 0")
            
            for rc in rcs:
                is_principal = (rc == principal)
                marker = " (PRINCIPAL)" if is_principal else " (non-principal)"
                
                print(f"\n    RC-graph{marker}:")
                print(f"      Weight: {rc.crystal_weight}")
                print(f"      Descents of perm: {rc.perm.descents()}")
                print(f"      i={i} in descents: {i in rc.perm.descents()}")
                
                # Check divdiff_desc - it returns a SET
                desc_result = rc.divdiff_desc(i)
                if len(desc_result) == 0:  # Empty set
                    print(f"      divdiff_desc({i}) = empty set (no root or not hw)")
                else:
                    print(f"      divdiff_desc({i}) = {len(desc_result)} element(s)")
                    for desc_rc in desc_result:
                        print(f"        perm {desc_rc.perm.trimcode}, weight {desc_rc.crystal_weight}")
                
                pretty_print(rc)
                print()
                
                # Show what dominant RCs this can extend with
                print(f"      Can extend products with dominant RCs where ∂_{i} ≠ 0:")
                extend_count = 0
                for dom_perm in dom_perms:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    if len(diff_dom) > 0:
                        extend_count += 1
                        if extend_count <= 3:  # Show first 3
                            print(f"        {dom_perm.trimcode}: ∂_{i} has {len(diff_dom)} terms")
                if extend_count > 3:
                    print(f"        ... and {extend_count - 3} more")
                print()
    
    print("\n" + "="*60)
    print("EXTENSION STRATEGY")
    print("="*60)
    print()
    print("For these pairs (rc, i), we can extend products via Leibniz:")
    print("  ∂_i(rc % dom) = s_i(rc) % ∂_i(dom)")
    print()
    print("Since ∂_i(rc) = 0, this simplifies to:")
    print("  ∂_i(rc % dom) = s_i(rc) % ∂_i(dom)")
    print()
    print("This allows us to define products rc % (non-dominant) by:")
    print("  rc % ∂_i(dom) := ∂_i(rc % dom) / s_i(rc)")
    print()
    print("where 'division' needs to be understood via ring localization.")
    print()
    
    return useful_pairs


def compare_principal_vs_all_rcs(n):
    """
    Compare: how many useful pairs for principal RCs vs all RCs?
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"PRINCIPAL vs ALL RC-GRAPHS (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant]
    rc_ring = RCGraphRing()
    
    # Count useful pairs for principal RCs only
    principal_useful = []
    for perm in perms:
        rc = RCGraph.principal_rc(perm, n-1)
        for i in range(1, n):
            diff_rc = rc_ring(rc).divdiff(i)
            if len(diff_rc) == 0:
                # Check if useful
                for dom_perm in dom_perms:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    if len(rc_ring(dom_rc).divdiff(i)) > 0:
                        principal_useful.append((rc, i))
                        break
    
    # Count useful pairs for all RCs using library function
    all_rcs = []
    for perm in perms:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    all_useful = []
    for rc in all_rcs:
        for i in range(1, n):
            diff_rc = rc_ring(rc).divdiff(i)
            if len(diff_rc) == 0:
                for dom_perm in dom_perms:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    if len(rc_ring(dom_rc).divdiff(i)) > 0:
                        all_useful.append((rc, i))
                        break
    
    print(f"\nStatistics:")
    print(f"  Principal RCs only: {len(principal_useful)} useful pairs")
    print(f"  All RCs:            {len(all_useful)} useful pairs")
    print(f"  Gain from non-principal: {len(all_useful) - len(principal_useful)} additional pairs")
    print()
    
    # Find ALL examples where non-principal RC has zero divdiff but principal doesn't
    print("="*60)
    print("ALL CASES: Non-principal RC has zero divdiff but principal doesn't")
    print("="*60)
    print()
    
    principal_divdiff = {}
    for perm in perms:
        principal = RCGraph.principal_rc(perm, n-1)
        principal_divdiff[perm] = {}
        for i in range(1, n):
            principal_divdiff[perm][i] = len(rc_ring(principal).divdiff(i))
    
    count = 0
    for rc, i in all_useful:
        principal = RCGraph.principal_rc(rc.perm, n-1)
        if rc != principal and principal_divdiff[rc.perm][i] > 0:
            count += 1
            print(f"\n{count}. Permutation {rc.perm.trimcode}, index i={i}:")
            print(f"   Descents: {rc.perm.descents()}")
            print()
            print(f"   Principal RC-graph: divdiff({i}) ≠ 0 ({principal_divdiff[rc.perm][i]} terms)")
            print(f"   Weight: {principal.crystal_weight}")
            pretty_print(principal)
            print()
            print(f"   This non-principal RC: divdiff({i}) = 0")
            print(f"   Weight: {rc.crystal_weight}")
            pretty_print(rc)
            print()
            print(f"   INSIGHT: Same permutation, different RC-graphs → different divdiff!")
            print(f"   This is because:")
            if i in rc.perm.descents():
                print(f"     - i={i} IS a descent")
                print(f"     - But root structure differs between principal and non-principal")
                print(f"     - Non-principal: α_i not in right place or deletion doesn't give hw")
            print()
    
    print(f"\nTotal: Found {count} cases where non-principal RC has zero divdiff")
    print(f"but principal RC does not.")
    print()
    print("This demonstrates that RC-graphs with the SAME permutation")
    print("can have DIFFERENT divided difference behavior depending on")
    print("their crystal structure!")
    print()


def analyze_s_i_action(n):
    """
    Analyze the s_i action on RC-graphs.
    
    For the Leibniz formula:
        ∂_i(a % b) = ∂_i(a) % b + s_i(a) % ∂_i(b)
    
    We need to understand what s_i(a) means.
    
    From Demazure operator theory:
        - If i ∈ Des(a.perm): s_i(a) acts via ∂_i
        - If i ∉ Des(a.perm): s_i(a) = left multiplication by s_i
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"UNDERSTANDING s_i ACTION (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    print("\nThe s_i action on RC-graphs (Demazure operators):")
    print()
    print("Case 1: i ∈ Des(a.perm)")
    print("  s_i(a) = ∂_i(a)  (Demazure relation)")
    print("  This is a sum of RC-graphs from the i-string")
    print()
    print("Case 2: i ∉ Des(a.perm)")
    print("  s_i(a) = s_i · a  (left multiplication)")
    print("  where s_i · a means the RC-graph for s_i * a.perm")
    print("  This should be a single RC-graph")
    print()
    
    # Test this on ALL permutations
    print("="*60)
    print("COMPLETE ANALYSIS OF s_i ACTION")
    print("="*60)
    print()
    
    for perm in perms:
        rc = RCGraph.principal_rc(perm, n-1)
        
        print(f"\nPermutation {perm.trimcode}:")
        print(f"  Descents: {perm.descents()}")
        print(f"  Weight: {rc.crystal_weight}")
        
        for i in range(1, n):
            s_i_perm = Permutation.ref_product(i)
            
            if i in perm.descents():
                # Case 1: s_i acts as ∂_i
                s_i_a = rc_ring(rc).divdiff(i)
                print(f"\n  i={i} (DESCENT):")
                print(f"    s_{i}(a) = ∂_{i}(a) = {len(s_i_a)} terms")
                if len(s_i_a) > 0:
                    print(f"    Terms:")
                    for term_rc, coeff in s_i_a.items():
                        print(f"      Coeff {coeff}, perm {term_rc.perm.trimcode}, weight {term_rc.crystal_weight}")
            else:
                # Case 2: s_i multiplies on left
                s_i_a_perm = s_i_perm * perm
                s_i_a_rc = RCGraph.principal_rc(s_i_a_perm, n-1)
                print(f"\n  i={i} (NOT descent):")
                print(f"    s_{i}(a) = s_{i} · a = single RC with perm {s_i_a_perm.trimcode}")
                print(f"    Weight: {s_i_a_rc.crystal_weight}")
        print()


def implement_leibniz_extension_correct(n):
    """
    Implement Leibniz extension with correct understanding.
    
    When divdiff(a, i) = 0, the Leibniz formula gives:
        ∂_i(a % dom) = s_i(a) % ∂_i(dom)
    
    Where:
    - ∂_i(dom) is a sum over the i-string of dom
    - s_i(a) is either ∂_i(a) (if descent) or s_i · a (if not)
    - The LHS is computable (we know a % dom)
    - The RHS defines products s_i(a) % (i-string elements)
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"IMPLEMENTING LEIBNIZ EXTENSION (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant]
    rc_ring = RCGraphRing()
    
    # Find ALL test cases including non-principal RCs
    test_cases = []
    
    # Get all RCs using library function
    all_rcs = []
    for perm in perms:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    for a_rc in all_rcs:
        for i in range(1, n):
            diff_a = rc_ring(a_rc).divdiff(i)
            
            if len(diff_a) == 0:
                # Found a useful case
                for dom_perm in dom_perms:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    
                    if len(diff_dom) > 0:
                        principal = RCGraph.principal_rc(a_rc.perm, n-1)
                        is_principal = (a_rc == principal)
                        test_cases.append((a_rc, dom_rc, i, is_principal))
                        break
                break
    
    print(f"\nFound {len(test_cases)} test cases (ALL RC-graphs with useful zero divdiff)")
    print()
    print("="*60)
    print("COMPLETE ANALYSIS OF ALL TEST CASES")
    print("="*60)
    
    for case_num, (a_rc, dom_rc, i, is_principal) in enumerate(test_cases, 1):
        marker = " (PRINCIPAL)" if is_principal else " (NON-PRINCIPAL)"
        print(f"\n\n{'='*60}")
        print(f"TEST CASE {case_num}{marker}")
        print(f"{'='*60}")
        print(f"\nLeft factor a: perm {a_rc.perm.trimcode}")
        print(f"Right factor dom: perm {dom_rc.perm.trimcode} (dominant)")
        print(f"Index: i={i}")
        print()
        
        print("Left factor RC-graph a:")
        print(f"  Weight: {a_rc.crystal_weight}")
        print(f"  Descents: {a_rc.perm.descents()}")
        print(f"  ∂_{i}(a) = 0 (verified)")
        pretty_print(a_rc)
        print()
        
        print("Right factor RC-graph dom:")
        print(f"  Weight: {dom_rc.crystal_weight}")
        pretty_print(dom_rc)
        print()
        
        # LHS: ∂_i(a % dom)
        print(f"Computing LHS: ∂_{i}(a % dom)")
        prod_a_dom = rc_ring(a_rc) % rc_ring(dom_rc)
        print(f"  a % dom has {len(prod_a_dom)} terms")
        lhs = prod_a_dom.divdiff(i)
        print(f"  ∂_{i}(a % dom) has {len(lhs)} terms")
        
        if len(lhs) > 0:
            print(f"\n  LHS terms:")
            for lhs_rc, lhs_coeff in lhs.items():
                print(f"    Coeff {lhs_coeff}, perm {lhs_rc.perm.trimcode}, weight {lhs_rc.crystal_weight}")
        print()
        
        # RHS: s_i(a) % ∂_i(dom)
        print(f"Computing RHS: s_{i}(a) % ∂_{i}(dom)")
        
        # Compute s_i(a)
        s_i_perm = Permutation.ref_product(i)
        
        if i in a_rc.perm.descents():
            s_i_a = rc_ring(a_rc).divdiff(i)
            print(f"  s_{i}(a) = ∂_{i}(a) has {len(s_i_a)} terms (descent case)")
            if len(s_i_a) > 0:
                print(f"  But wait - we said ∂_{i}(a) = 0!")
                print(f"  This should not happen!")
        else:
            s_i_a_perm = s_i_perm * a_rc.perm
            s_i_a_rc = RCGraph.principal_rc(s_i_a_perm, n-1)
            s_i_a = rc_ring(s_i_a_rc)
            print(f"  s_{i}(a) = s_{i} · a = single RC with perm {s_i_a_perm.trimcode}")
            print(f"  Weight: {s_i_a_rc.crystal_weight}")
            print(f"  (non-descent case)")
        
        # Compute ∂_i(dom)
        diff_dom = rc_ring(dom_rc).divdiff(i)
        print(f"\n  ∂_{i}(dom) has {len(diff_dom)} terms (i-string)")
        print(f"  ∂_{i}(dom) terms:")
        for dom_i_rc, dom_i_coeff in diff_dom.items():
            print(f"    Coeff {dom_i_coeff}, perm {dom_i_rc.perm.trimcode}, weight {dom_i_rc.crystal_weight}")
        
        # RHS should be s_i(a) % ∂_i(dom)
        expected_products = len(s_i_a) * len(diff_dom)
        
        print(f"\n  RHS involves {len(s_i_a)} × {len(diff_dom)} = {expected_products} products")
        print(f"  These products are what we want to DEFINE using the localization!")
        print()
        
        # Match by crystal weights
        lhs_weights = {}
        for rc, coeff in lhs.items():
            weight = tuple(rc.crystal_weight)
            if weight not in lhs_weights:
                lhs_weights[weight] = []
            lhs_weights[weight].append((rc, coeff))
        
        print(f"  LHS has {len(lhs_weights)} distinct weights:")
        for weight, terms in lhs_weights.items():
            print(f"    Weight {weight}: {len(terms)} term(s)")
        print()
        
        print("  LOCALIZATION STRATEGY:")
        print(f"    We know: LHS = ∂_{i}(a % dom) = s_{i}(a) % ∂_{i}(dom) = RHS")
        print(f"    We know: a % dom (computed above)")
        print(f"    We know: s_{i}(a) (computed above)")
        print(f"    We know: ∂_{i}(dom) (computed above)")
        print(f"    We want: To define products s_{i}(a) % each term of ∂_{i}(dom)")
        print(f"    Method: 'Divide' LHS by s_{i}(a) to extract the products with ∂_{i}(dom) terms")
        print()


def explore_ring_localization(n):
    """
    Explore the ring-theoretic structure and potential localization.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"RING LOCALIZATION AND EXTENSION (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant]
    rc_ring = RCGraphRing()
    
    print("\n1. STRUCTURE OF SCHUBERT POLYNOMIAL RING")
    print("="*60)
    print(f"   {len(dom_perms)} dominant permutations = Schubert basis")
    print(f"   Dominant perms: {[p.trimcode for p in dom_perms]}")
    print("   This forms a commutative ring (via % product)")
    print("   No zero divisors (weight is additive)")
    print()
    
    # Test: find kernel of divided differences on principal RCs
    print("2. KERNEL OF DIVIDED DIFFERENCES")
    print("="*60)
    
    for i in range(1, n):
        ker_i = []
        for perm in perms:
            rc = RCGraph.principal_rc(perm, n-1)
            diff = rc_ring(rc).divdiff(i)
            if len(diff) == 0:
                ker_i.append(perm)
        
        print(f"\n   ker(∂_{i}) on principal RCs: {len(ker_i)} elements")
        print(f"   Elements: {[p.trimcode for p in ker_i]}")
        
        # These are exactly the perms with no descent at i
        non_descents = [p for p in perms if i not in p.descents()]
        print(f"   Non-descents at i: {len(non_descents)} elements")
        print(f"   Elements: {[p.trimcode for p in non_descents]}")
        print(f"   Match? {set(ker_i) == set(non_descents)}")
    
    print("\n3. WEIGHT-PRESERVING EXTENSION")
    print("="*60)
    print("   Goal: Define products a % b where b is non-dominant")
    print("   such that weight(a % b) = weight(a) + weight(b)")
    print()
    
    # Get all RC-graphs
    all_rcs = []
    for perm in perms:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    print(f"   Total RC-graphs: {len(all_rcs)}")
    print(f"   Principal RCs: {len(perms)}")
    print(f"   Non-principal: {len(all_rcs) - len(perms)}")
    print()
    
    # For ALL non-principal RCs, show crystal weight
    print("   ALL non-principal RC-graphs and their weights:")
    print()
    non_principal_count = 0
    for rc in all_rcs:
        principal = RCGraph.principal_rc(rc.perm, n-1)
        if rc != principal:
            non_principal_count += 1
            print(f"\n   {non_principal_count}. Perm {rc.perm.trimcode}:")
            print(f"      Weight: {rc.crystal_weight}")
            print(f"      Descents: {rc.perm.descents()}")
            pretty_print(rc)
    
    print(f"\n\n   Summary: {non_principal_count} non-principal RC-graphs found")
    print()


def compute_new_products_via_leibniz(n):
    """
    Actually compute NEW products that aren't (any RC) % (dominant RC).
    
    Strategy:
    1. Find (a, i) where ∂_i(a) = 0
    2. Find dominant dom where ∂_i(dom) ≠ 0
    3. Compute LHS = ∂_i(a % dom)
    4. Know RHS = s_i(a) % ∂_i(dom)
    5. Each term in ∂_i(dom) should be NON-DOMINANT
    6. Use this to deduce: s_i(a) % (non-dominant term)
    
    EXHAUSTIVE SEARCH: Try all possible opportunities.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"COMPUTING NEW PRODUCTS VIA LEIBNIZ (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant]
    rc_ring = RCGraphRing()
    
    # Get all RCs
    all_rcs = []
    for perm in perms:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    print(f"\nTotal permutations: {len(perms)}")
    print(f"Dominant permutations: {len(dom_perms)}")
    print(f"  Dominant: {[p.trimcode for p in dom_perms]}")
    print(f"Total RC-graphs: {len(all_rcs)}")
    print()
    
    new_products_computed = []
    opportunities_tried = 0
    
    print("Exhaustive search for opportunities to compute new products...")
    print()
    
    for a_rc in all_rcs:
        for i in range(1, n):
            diff_a = rc_ring(a_rc).divdiff(i)
            
            if len(diff_a) == 0:  # Found useful a
                # Find dominant where ∂_i ≠ 0
                for dom_perm in dom_perms:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    
                    if len(diff_dom) > 0:
                        opportunities_tried += 1
                        
                        # Check if ANY term in ∂_i(dom) is truly non-dominant
                        has_non_dominant = False
                        for b_rc, b_coeff in diff_dom.items():
                            if not b_rc.perm.is_dominant:
                                has_non_dominant = True
                                break
                        
                        if not has_non_dominant:
                            # All terms are dominant, skip
                            continue
                        
                        # We have a real opportunity!
                        print("="*60)
                        print(f"OPPORTUNITY {opportunities_tried}")
                        print("="*60)
                        print(f"\nLeft factor a: perm {a_rc.perm.trimcode}, weight {a_rc.crystal_weight}")
                        principal_a = RCGraph.principal_rc(a_rc.perm, n-1)
                        is_principal_a = (a_rc == principal_a)
                        print(f"  Is principal? {is_principal_a}")
                        pretty_print(a_rc)
                        print()
                        
                        print(f"Right factor dom: perm {dom_rc.perm.trimcode}, weight {dom_rc.crystal_weight}")
                        print(f"  Is dominant? {dom_rc.perm.is_dominant}")
                        pretty_print(dom_rc)
                        print()
                        
                        print(f"Index i={i}")
                        print(f"  ∂_{i}(a) = 0 (verified)")
                        print(f"  i in Des(a)? {i in a_rc.perm.descents()}")
                        print()
                        
                        # Compute known LHS
                        prod_a_dom = rc_ring(a_rc) % rc_ring(dom_rc)
                        lhs = prod_a_dom.divdiff(i)
                        
                        print(f"Known product: a % dom")
                        print(f"  Has {len(prod_a_dom)} terms:")
                        for prod_rc, prod_coeff in prod_a_dom.items():
                            print(f"    Coeff {prod_coeff}: perm {prod_rc.perm.trimcode}, weight {prod_rc.crystal_weight}")
                        print()
                        
                        print(f"Known: ∂_{i}(a % dom) = LHS has {len(lhs)} terms:")
                        for lhs_rc, lhs_coeff in lhs.items():
                            print(f"  Coeff {lhs_coeff}: perm {lhs_rc.perm.trimcode}, weight {lhs_rc.crystal_weight}")
                        print()
                        
                        # Compute s_i(a)
                        if i in a_rc.perm.descents():
                            print(f"ERROR: i={i} is a descent of a, but ∂_i(a) = 0")
                            print(f"This should not happen with principal RCs!")
                            print(f"But can happen with non-principal RCs.")
                            print()
                            continue
                        
                        s_i_perm = Permutation.ref_product(i)
                        s_i_a_perm = s_i_perm * a_rc.perm
                        s_i_a_rc = RCGraph.principal_rc(s_i_a_perm, n-1)
                        
                        print(f"s_{i}(a) = s_{i} · a:")
                        print(f"  Single RC with perm {s_i_a_perm.trimcode}, weight {s_i_a_rc.crystal_weight}")
                        pretty_print(s_i_a_rc)
                        print()
                        
                        # The i-string ∂_i(dom)
                        print(f"∂_{i}(dom) has {len(diff_dom)} terms:")
                        for b_rc, b_coeff in diff_dom.items():
                            print(f"  Coeff {b_coeff}: perm {b_rc.perm.trimcode}, weight {b_rc.crystal_weight}")
                            print(f"    Is dominant? {b_rc.perm.is_dominant}")
                        print()
                        
                        # Now: LHS = s_i(a) % ∂_i(dom)
                        print("LEIBNIZ EQUATION:")
                        print(f"  ∂_{i}(a % dom) = s_{i}(a) % ∂_{i}(dom)")
                        print()
                        
                        # Try to deduce the new products
                        print("ATTEMPTING TO DEDUCE NEW PRODUCTS:")
                        print()
                        
                        # For each term b in ∂_i(dom), we want: s_i(a) % b
                        for b_rc, b_coeff in diff_dom.items():
                            if b_rc.perm.is_dominant:
                                print(f"SKIPPING: {b_rc.perm.trimcode} is DOMINANT")
                                print(f"  This product is already known")
                                print()
                                continue
                            
                            print(f"NEW PRODUCT: s_{i}(a) % {b_rc.perm.trimcode}")
                            print(f"  Left:  perm {s_i_a_perm.trimcode}, weight {s_i_a_rc.crystal_weight}")
                            print(f"  Right: perm {b_rc.perm.trimcode}, weight {b_rc.crystal_weight} (NON-DOMINANT!)")
                            print(f"  Right descents: {b_rc.perm.descents()}")
                            print()
                            
                            # Expected weight
                            expected_weight = tuple([s_i_a_rc.crystal_weight[j] + b_rc.crystal_weight[j] 
                                                    for j in range(len(b_rc.crystal_weight))])
                            print(f"  Expected total weight: {expected_weight}")
                            print()
                            
                            # Check if we can compute this directly first
                            can_compute_directly = False
                            try:
                                direct_product = rc_ring(s_i_a_rc) % rc_ring(b_rc)
                                can_compute_directly = True
                                print(f"  DIRECT COMPUTATION WORKS! Result has {len(direct_product)} terms:")
                                for prod_rc, prod_coeff in direct_product.items():
                                    print(f"    Coeff {prod_coeff}: perm {prod_rc.perm.trimcode}, weight {prod_rc.crystal_weight}")
                                print()
                                
                                # Verify weight preservation
                                all_weights_match = True
                                for prod_rc, prod_coeff in direct_product.items():
                                    expected_w = [s_i_a_rc.crystal_weight[j] + b_rc.crystal_weight[j] 
                                                      for j in range(len(b_rc.crystal_weight))]
                                    actual_w = list(prod_rc.crystal_weight)
                                    if expected_w == actual_w:
                                        print(f"    ✓ Weight preserved for {prod_rc.perm.trimcode}")
                                    else:
                                        print(f"    ✗ Weight NOT preserved for {prod_rc.perm.trimcode}")
                                        print(f"      Expected: {expected_w}")
                                        print(f"      Actual:   {actual_w}")
                                        all_weights_match = False
                                
                                if all_weights_match:
                                    print("  ✓ All weights match - product is already correctly defined!")
                                print()
                                
                            except Exception as e:
                                print(f"  Direct computation FAILED: {e}")
                                print(f"  THIS IS A TRULY NEW PRODUCT WE NEED TO DEFINE!")
                                print()
                                
                                # Try to extract from LHS using weight matching
                                print(f"  Searching LHS for terms with weight {expected_weight}...")
                                matching_terms = []
                                for lhs_rc, lhs_coeff in lhs.items():
                                    if tuple(lhs_rc.crystal_weight) == expected_weight:
                                        matching_terms.append((lhs_rc, lhs_coeff))
                                        print(f"    Found: coeff {lhs_coeff}, perm {lhs_rc.perm.trimcode}")
                                
                                if matching_terms:
                                    print()
                                    print(f"  DEDUCED PRODUCT (divided by b_coeff={b_coeff}):")
                                    for lhs_rc, lhs_coeff in matching_terms:
                                        deduced_coeff = lhs_coeff / b_coeff
                                        print(f"    Coeff {deduced_coeff}: perm {lhs_rc.perm.trimcode}, weight {lhs_rc.crystal_weight}")
                                    
                                    new_products_computed.append({
                                        'left': (s_i_a_perm.trimcode, tuple(s_i_a_rc.crystal_weight)),
                                        'right': (b_rc.perm.trimcode, tuple(b_rc.crystal_weight)),
                                        'result': [(rc.perm.trimcode, tuple(rc.crystal_weight), coeff) 
                                                  for rc, coeff in matching_terms],
                                        'method': 'leibniz_localization',
                                        'can_verify_directly': False
                                    })
                                else:
                                    print(f"  ✗ No matching terms found in LHS!")
                                print()
                            
                            if can_compute_directly:
                                # Still record that this product is computable
                                print("  This product is already computable, not truly new.")
                                print()
    
    print("\n" + "="*60)
    print("SUMMARY OF NEW PRODUCTS COMPUTED")
    print("="*60)
    print()
    print(f"Total opportunities examined: {opportunities_tried}")
    print(f"Total NEW products deduced (not directly computable): {len(new_products_computed)}")
    print()
    
    if len(new_products_computed) == 0:
        print("No truly new products found!")
        print("This means the product is already fully defined for all RC-graph pairs,")
        print("or we need to look at higher n to find undefined products.")
    else:
        for idx, prod_info in enumerate(new_products_computed, 1):
            print(f"{idx}. {prod_info['left'][0]} % {prod_info['right'][0]}")
            print(f"   Left weight:  {prod_info['left'][1]}")
            print(f"   Right weight: {prod_info['right'][1]} (NON-DOMINANT)")
            print(f"   Result: {len(prod_info['result'])} term(s)")
            for perm_str, weight, coeff in prod_info['result']:
                print(f"     Coeff {coeff}: perm {perm_str}, weight {weight}")
            print()
    
    return new_products_computed


if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    
    clarify_divdiff_structure()
    zero_cases = when_is_divdiff_zero(n)
    compare_principal_vs_all_rcs(n)
    useful_pairs = find_useful_zero_divdiff_pairs(n)
    analyze_s_i_action(n)
    implement_leibniz_extension_correct(n)
    explore_ring_localization(n)
    
    # NEW: Actually compute new products!
    new_products = compute_new_products_via_leibniz(n)
    
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    print()
    print(f"Total permutations (n={n}): {len(Permutation.all_permutations(n))}")
    print(f"Total RC-graphs: {sum(len(RCGraph.all_rc_graphs(p, n-1)) for p in Permutation.all_permutations(n))}")
    print(f"Total useful (rc, i) pairs for extension: {len(useful_pairs)}")
    print(f"NEW PRODUCTS COMPUTED (truly new, not directly computable): {len(new_products)}")
    print()
    print("Key insights:")
    print("  1. Different RC-graphs with SAME permutation have DIFFERENT divdiff")
    print("  2. Non-principal RCs give additional useful pairs for Leibniz extension")
    print("  3. Need S_∞ action on non-principal RCs for complete theory")
    print("  4. Localization approach allows systematic product extension")
    if len(new_products) > 0:
        print(f"  5. Successfully deduced {len(new_products)} new products via Leibniz!")
    else:
        print(f"  5. All products appear to be already defined (try higher n)")