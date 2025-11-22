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
        'not_hw': [],
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
                    # i is a descent, check if α_i is in row i
                    desc_result = rc.divdiff_desc(i)
                    
                    if desc_result is None:
                        # α_i not in row i, or doesn't yield hw
                        zero_cases['no_root'].append((rc, i))
                    else:
                        # Shouldn't happen if divdiff = 0 but divdiff_desc ≠ 0
                        print(f"  ⚠ Unexpected: RC-graph with perm {rc.perm.trimcode}, i={i}")
                        print(f"    divdiff = 0 but divdiff_desc ≠ 0")
    
    print(f"\nCases where divdiff(rc, i) = 0:")
    print()
    print(f"1. No descent at i: {len(zero_cases['no_descent'])} cases")
    print(f"2. Descent but no α_i in row i: {len(zero_cases['no_root'])} cases")
    print()
    
    # Show examples
    print("Examples of 'no descent' cases:")
    for rc, i in zero_cases['no_descent'][:3]:
        print(f"  RC {rc.perm.trimcode}: i={i} (not in {rc.perm.descents()})")
        pretty_print(rc)
        print()
    
    print("\nExamples of 'descent but no root' cases:")
    for rc, i in zero_cases['no_root'][:3]:
        print(f"  RC {rc.perm.trimcode}: i={i} (in {rc.perm.descents()} but α_i not in row i)")
        pretty_print(rc)
        print()
    
    return zero_cases


def find_useful_zero_divdiff_pairs(n):
    """
    Find pairs (rc, i) where divdiff(rc, i) = 0 AND this is useful
    for extending products via Leibniz rule.
    
    KEY: Check ALL RC-graphs, not just principal ones!
    
    These are cases where we can use:
        ∂_i(rc % dom) = s_i(rc) % ∂_i(dom)
    
    to extend the product to non-dominant right factors.
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
    
    print(f"\nFound {len(useful_pairs)} useful pairs where:")
    print("  - divdiff(rc, i) = 0")
    print("  - ∃ dominant dom with divdiff(dom, i) ≠ 0")
    print()
    
    # Group by permutation
    by_perm = {}
    for rc, i in useful_pairs:
        if rc.perm not in by_perm:
            by_perm[rc.perm] = {}
        if i not in by_perm[rc.perm]:
            by_perm[rc.perm][i] = []
        by_perm[rc.perm][i].append(rc)
    
    print("Useful pairs by permutation and index (showing first 5 perms):")
    for perm in sorted(by_perm.keys(), key=lambda p: p.inv)[:5]:
        print(f"\n  Permutation {perm.trimcode}:")
        for i in sorted(by_perm[perm].keys()):
            rcs = by_perm[perm][i]
            principal = RCGraph.principal_rc(perm, n-1)
            if len(rcs) == 1 and rcs[0] == principal:
                print(f"    i={i}: principal RC only")
                pretty_print(principal)
                print()
            else:
                print(f"    i={i}: {len(rcs)} RC-graphs")
                for rc in rcs[:2]:  # Show first 2
                    is_principal = (rc == principal)
                    marker = " (principal)" if is_principal else ""
                    print(f"      - weight={rc.crystal_weight}{marker}")
                    pretty_print(rc)
                    print()
    
    print("\nFor these pairs, we can extend:")
    print("  rc % ∂_i(dom)  is determined by  ∂_i(rc % dom)")
    print("  via the Leibniz formula:")
    print("    ∂_i(rc % dom) = s_i(rc) % ∂_i(dom)")
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
    
    print(f"\nPrincipal RCs only: {len(principal_useful)} useful pairs")
    print(f"All RCs:            {len(all_useful)} useful pairs")
    print(f"Gain from non-principal: {len(all_useful) - len(principal_useful)} additional pairs")
    print()
    
    # Find examples where non-principal RC has zero divdiff but principal doesn't
    print("Examples where non-principal RC gives zero but principal doesn't:")
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
            if count <= 3:
                print(f"  Perm {rc.perm.trimcode}, i={i}:")
                print(f"    Principal: divdiff ≠ 0 ({principal_divdiff[rc.perm][i]} terms)")
                pretty_print(principal)
                print(f"    This RC:   divdiff = 0")
                pretty_print(rc)
                print()
    
    print(f"\nFound {count} cases where non-principal RC has zero divdiff")
    print(f"but principal RC does not.")
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
    
    print("\nThe s_i action on RC-graphs:")
    print()
    print("Case 1: i ∈ Des(a.perm)")
    print("  s_i(a) = ∂_i(a)  (Demazure relation)")
    print()
    print("Case 2: i ∉ Des(a.perm)")
    print("  s_i(a) = s_i · a  (left multiplication)")
    print("  where s_i · a means the RC-graph for s_i * a.perm")
    print()
    
    # Test this on examples
    print("Examples:")
    print()
    
    for perm in perms[:3]:
        rc = RCGraph.principal_rc(perm, n-1)
        
        for i in range(1, n):
            s_i_perm = Permutation.ref_product(i)
            
            if i in perm.descents():
                # Case 1: s_i acts as ∂_i
                s_i_a = rc_ring(rc).divdiff(i)
                print(f"  {perm.trimcode}, i={i} (descent):")
                print(f"    s_{i}(a) = ∂_{i}(a) = {len(s_i_a)} terms")
            else:
                # Case 2: s_i multiplies on left
                s_i_a_perm = s_i_perm * perm
                s_i_a_rc = RCGraph.principal_rc(s_i_a_perm, n-1)
                print(f"  {perm.trimcode}, i={i} (not descent):")
                print(f"    s_{i}(a) = s_{i} · a = {s_i_a_perm.trimcode}")


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
    
    # Find test cases including non-principal RCs
    test_cases = []
    
    # Get all RCs using library function (limit to first 10 perms for speed)
    all_rcs = []
    for perm in perms[:10]:
        all_rcs.extend(RCGraph.all_rc_graphs(perm, n-1))
    
    for a_rc in all_rcs[:20]:
        for i in range(1, n):
            diff_a = rc_ring(a_rc).divdiff(i)
            
            if len(diff_a) == 0:
                # Found a useful case
                for dom_perm in dom_perms[:2]:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    
                    if len(diff_dom) > 0:
                        principal = RCGraph.principal_rc(a_rc.perm, n-1)
                        is_principal = (a_rc == principal)
                        test_cases.append((a_rc, dom_rc, i, is_principal))
                        break
                break
        
        if len(test_cases) >= 3:
            break
    
    print(f"\nTesting {len(test_cases)} cases (including non-principal RCs):")
    print()
    
    for a_rc, dom_rc, i, is_principal in test_cases:
        marker = " (principal)" if is_principal else " (NON-PRINCIPAL)"
        print(f"Test: a={a_rc.perm.trimcode}{marker}, dom={dom_rc.perm.trimcode}, i={i}")
        print("RC-graph a:")
        pretty_print(a_rc)
        print()
        
        # LHS: ∂_i(a % dom)
        prod_a_dom = rc_ring(a_rc) % rc_ring(dom_rc)
        lhs = prod_a_dom.divdiff(i)
        
        print(f"  LHS = ∂_{i}(a % dom) has {len(lhs)} terms")
        
        # RHS: s_i(a) % ∂_i(dom)
        
        # Compute s_i(a)
        s_i_perm = Permutation.ref_product(i)
        
        if i in a_rc.perm.descents():
            s_i_a = rc_ring(a_rc).divdiff(i)
            print(f"  s_{i}(a) = ∂_{i}(a) has {len(s_i_a)} terms (descent case)")
        else:
            s_i_a_perm = s_i_perm * a_rc.perm
            s_i_a_rc = RCGraph.principal_rc(s_i_a_perm, n-1)
            s_i_a = rc_ring(s_i_a_rc)
            print(f"  s_{i}(a) = s_{i} · a = {s_i_a_perm.trimcode} (non-descent case)")
        
        # Compute ∂_i(dom)
        diff_dom = rc_ring(dom_rc).divdiff(i)
        print(f"  ∂_{i}(dom) has {len(diff_dom)} terms (i-string)")
        
        # RHS should be s_i(a) % ∂_i(dom)
        # This is a sum of products
        expected_products = len(s_i_a) * len(diff_dom)
        
        print(f"  RHS should involve {expected_products} products")
        print(f"  These products are what we want to DEFINE")
        print()
        
        # Match by crystal weights
        lhs_weights = {}
        for rc, coeff in lhs.items():
            weight = tuple(rc.crystal_weight)
            if weight not in lhs_weights:
                lhs_weights[weight] = []
            lhs_weights[weight].append((rc, coeff))
        
        print(f"  LHS has {len(lhs_weights)} distinct weights")
        print()


if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    
    clarify_divdiff_structure()
    
    zero_cases = when_is_divdiff_zero(n)
    
    if n <= 4:
        compare_principal_vs_all_rcs(n)
    
    useful_pairs = find_useful_zero_divdiff_pairs(n)
    
    analyze_s_i_action(n)
    
    if n <= 4:
        implement_leibniz_extension_correct(n)
    
    print("\n" + "="*80)
    print("SUMMARY: CORRECTED UNDERSTANDING")
    print("="*80)
    print()
    print("divdiff_desc(rc, i) ≠ 0 requires:")
    print("  1. i ∈ Des(rc.perm)")
    print("  2. α_i in row i of rc")
    print("  3. Deletion yields i-highest weight")
    print()
    print("divdiff(rc, i) = full i-string starting from divdiff_desc(rc, i)")
    print()
    print("IMPORTANT: Different RC-graphs with the same permutation")
    print("           can have different divdiff behavior!")
    print()
    print("For Leibniz extension when divdiff(rc, i) = 0:")
    print("  ∂_i(rc % dom) = s_i(rc) % ∂_i(dom)")
    print()
    print("Where s_i(rc) is:")
    print("  - ∂_i(rc) if i is a descent")
    print("  - s_i · rc otherwise")
    print()
    
    # Display all useful pairs
    print("\n" + "="*80)
    print("ALL USEFUL PAIRS FOR PRODUCT EXTENSION")
    print("="*80)
    print()
    print(f"Total: {len(useful_pairs)} pairs (rc, i) where:")
    print("  - divdiff(rc, i) = 0")
    print("  - Can extend product to non-dominant right factors")
    print()
    
    # Group by permutation and index
    by_perm = {}
    for rc, i in useful_pairs:
        if rc.perm not in by_perm:
            by_perm[rc.perm] = {}
        if i not in by_perm[rc.perm]:
            by_perm[rc.perm][i] = []
        by_perm[rc.perm][i].append(rc)
    
    print("Complete list (grouped by permutation, showing first 5):")
    for perm in sorted(by_perm.keys(), key=lambda p: (p.inv, p.trimcode))[:5]:
        print(f"\n  Permutation {perm.trimcode}:")
        principal = RCGraph.principal_rc(perm, n-1)
        for i in sorted(by_perm[perm].keys()):
            rcs = by_perm[perm][i]
            
            if len(rcs) == 1 and rcs[0] == principal:
                print(f"    i={i}: principal RC only")
                pretty_print(principal)
                print()
            else:
                print(f"    i={i}: {len(rcs)} RC-graphs (showing first 2):")
                for rc in rcs[:2]:
                    is_principal = (rc == principal)
                    marker = " (principal)" if is_principal else ""
                    print(f"      weight={rc.crystal_weight}{marker}")
                    pretty_print(rc)
                    print()
    
    print()
    print("These pairs allow extending the % product via:")
    print("  rc % ∂_i(dom) := ∂_i(rc % dom) / s_i(rc)")
    print()
    
    # Count non-principal useful pairs
    non_principal_count = sum(1 for rc, i in useful_pairs 
                              if rc != RCGraph.principal_rc(rc.perm, n-1))
    print(f"By checking ALL RC-graphs (not just principal),")
    print(f"we found {non_principal_count} additional pairs!")