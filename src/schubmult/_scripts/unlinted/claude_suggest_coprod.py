from sympy import pretty_print

def test_local_canonicality_conditions(n):
    """
    Test various LOCAL conditions on RC-graphs that might ensure:
    1. Canonical choice within each coproduct
    2. Automatic disjointness between different coproducts
    
    Key insight: The condition should depend only on properties of the
    individual RC-graph pair (rc1, rc2), not on what other coproducts exist.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"TESTING LOCAL CANONICALITY CONDITIONS (n={n})")
    print("="*80)
    
    print("\nGoal: Find a LOCAL condition that:")
    print("  1. Uniquely selects one (rc1, rc2) pair for each (perm1, perm2)")
    print("  2. Depends only on properties of (rc1, rc2) itself")
    print("  3. Automatically ensures disjoint coproducts")
    print()
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    # Collect ALL coproduct terms from ALL permutations
    all_coprod_terms = {}  # Maps (rc1, rc2) -> list of source perms
    
    for perm in perms:
        principal = RCGraph.principal_rc(perm, n-1)
        coprod = rc_ring.coproduct_on_basis(principal)
        
        for (rc1, rc2), coeff in coprod.items():
            key = (rc1, rc2)
            if key not in all_coprod_terms:
                all_coprod_terms[key] = []
            all_coprod_terms[key].append(perm)
    
    print(f"Total unique (rc1, rc2) pairs across all coproducts: {len(all_coprod_terms)}")
    print()
    
    # Find overlapping terms
    overlapping_terms = {k: v for k, v in all_coprod_terms.items() if len(v) > 1}
    
    print(f"RC-graph pairs appearing in MULTIPLE coproducts: {len(overlapping_terms)}")
    
    if len(overlapping_terms) == 0:
        print("✓ ALREADY DISJOINT! The current coproduct is perfect!")
        return
    
    print(f"⚠ Found {len(overlapping_terms)} overlapping pairs")
    print()
    
    # Analyze the overlapping terms
    print("="*60)
    print("ANALYZING OVERLAPPING TERMS")
    print("="*60)
    
    for idx, ((rc1, rc2), source_perms) in enumerate(list(overlapping_terms.items())[:5]):
        print(f"\nOverlap {idx+1}: RC pair appears in {len(source_perms)} coproducts")
        print(f"  Permutation pair: ({rc1.perm.trimcode}, {rc2.perm.trimcode})")
        print(f"  Source perms: {[p.trimcode for p in source_perms]}")
        print(f"  Left RC:  weight={rc1.crystal_weight}, len_vec={rc1.length_vector}")
        print(f"  Right RC: weight={rc2.crystal_weight}, len_vec={rc2.length_vector}")
        
        # Test LOCAL conditions that might distinguish which source it belongs to
        print(f"\n  Testing local conditions:")
        
        # Condition 1: Weight sum matches one source uniquely
        weight_sum = sum(rc1.crystal_weight) + sum(rc2.crystal_weight)
        print(f"    Total weight: {weight_sum}")
        for src_perm in source_perms:
            src_principal = RCGraph.principal_rc(src_perm, n-1)
            src_weight_sum = sum(src_principal.crystal_weight)
            print(f"      Source {src_perm.trimcode} has weight sum: {src_weight_sum}")
            if src_weight_sum == weight_sum:
                print(f"        ✓ MATCH!")
        
        # Condition 2: Length vector sum
        len_vec_sum = sum(rc1.length_vector) + sum(rc2.length_vector)
        print(f"    Total length: {len_vec_sum}")
        for src_perm in source_perms:
            src_principal = RCGraph.principal_rc(src_perm, n-1)
            src_len_sum = sum(src_principal.length_vector)
            print(f"      Source {src_perm.trimcode} has length sum: {src_len_sum}")
            if src_len_sum == len_vec_sum:
                print(f"        ✓ MATCH!")
        
        # Condition 3: Product of vertical cuts should reconstruct source
        print(f"    Vertical cut compatibility:")
        for src_perm in source_perms:
            src_principal = RCGraph.principal_rc(src_perm, n-1)
            
            # Check if rc1 * rc2 "grows to" src_principal
            # via the vertical cut structure
            compatible = check_vertical_cut_origin(rc1, rc2, src_principal)
            print(f"      Compatible with {src_perm.trimcode}? {compatible}")
    
    print()


def check_vertical_cut_origin(rc1, rc2, target):
    """
    Check if (rc1, rc2) could have originated from target's coproduct
    based on vertical cut structure.
    
    This is a LOCAL condition that depends only on rc1, rc2, and target.
    """
    # If lengths don't add up correctly, can't be from this target
    if len(rc1) + len(rc2) != len(target):
        return False
    
    # Check if weight vectors are compatible
    if len(rc1.crystal_weight) + len(rc2.crystal_weight) != len(target.crystal_weight):
        return False
    
    # More sophisticated check: can we reconstruct target's structure?
    # This would need to check the vertical cut decomposition
    # For now, just check basic compatibility
    
    return True


def propose_local_canonical_condition(n):
    """
    Propose a LOCAL condition that ensures canonicity and disjointness.
    
    Key idea: The condition should be intrinsic to (rc1, rc2) and should
    "know" which source permutation it came from.
    """
    print("\n" + "="*80)
    print(f"PROPOSING LOCAL CANONICAL CONDITION (n={n})")
    print("="*80)
    
    print("\nPROPOSALS FOR LOCAL CONDITIONS:")
    print()
    
    print("PROPOSAL A: Highest Weight Dominance")
    print("  For (rc1, rc2) to appear in coproduct of w:")
    print("  - Let hw1, hw2 be highest weights of rc1, rc2")
    print("  - Let hw_w be highest weight of principal RC for w")
    print("  - CONDITION: hw_w should be 'reachable' from (hw1, hw2)")
    print("    via crystal operations that preserve structure")
    print()
    
    print("PROPOSAL B: Vertical Cut Signature")
    print("  Each (rc1, rc2) has a unique 'signature' based on:")
    print("  - Sequence of vertical cuts in rc1 and rc2")
    print("  - This signature should uniquely determine which w it came from")
    print("  - CONDITION: Only include (rc1, rc2) if its signature matches")
    print("    the signature pattern of w")
    print()
    
    print("PROPOSAL C: Inversion Number Accounting")
    print("  For (rc1, rc2) to appear in coproduct of w:")
    print("  - Let inv1, inv2 be inversions in rc1.perm, rc2.perm")
    print("  - Let inv_w be inversions in w")
    print("  - CONDITION: inv1 + inv2 = inv_w")
    print("  - PLUS: Additional structure on HOW inversions distribute")
    print()
    
    print("PROPOSAL D: Maximal Element Condition")
    print("  For each (perm1, perm2) pair, choose (rc1, rc2) that is:")
    print("  - Maximal in crystal order")
    print("  - Highest weight in the crystal component")
    print("  - This is LOCAL: depends only on (rc1, rc2), not other terms")
    print()
    
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    print("="*60)
    print("TESTING PROPOSAL C: INVERSION ACCOUNTING")
    print("="*60)
    print()
    
    # Test Proposal C
    all_valid_by_inv = {}  # Maps (rc1, rc2) -> list of valid sources by inversion rule
    
    for perm in perms:
        principal = RCGraph.principal_rc(perm, n-1)
        coprod = rc_ring.coproduct_on_basis(principal)
        inv_w = perm.inv
        
        for (rc1, rc2), coeff in coprod.items():
            inv1 = rc1.perm.inv
            inv2 = rc2.perm.inv
            
            key = (rc1, rc2)
            if key not in all_valid_by_inv:
                all_valid_by_inv[key] = []
            
            # Check if inversion rule is satisfied
            if inv1 + inv2 == inv_w:
                all_valid_by_inv[key].append(('inv_match', perm))
            else:
                all_valid_by_inv[key].append(('inv_mismatch', perm))
    
    # Count how many satisfy inversion rule
    inv_unique = sum(1 for terms in all_valid_by_inv.values() 
                     if len([t for t in terms if t[0] == 'inv_match']) == 1)
    
    print(f"Pairs with UNIQUE inversion match: {inv_unique} / {len(all_valid_by_inv)}")
    print()
    
    # Show examples where inversion rule works
    inv_works = {k: v for k, v in all_valid_by_inv.items() 
                 if len([t for t in v if t[0] == 'inv_match']) == 1}
    
    if inv_works:
        print(f"✓ Inversion rule uniquely identifies {len(inv_works)} pairs!")
        print(f"  Examples (first 3):")
        for idx, ((rc1, rc2), terms) in enumerate(list(inv_works.items())[:3]):
            matching = [t[1] for t in terms if t[0] == 'inv_match'][0]
            print(f"    Pair ({rc1.perm.trimcode}, {rc2.perm.trimcode})")
            print(f"      Inversions: {rc1.perm.inv} + {rc2.perm.inv} = {rc1.perm.inv + rc2.perm.inv}")
            print(f"      Uniquely matches source: {matching.trimcode} (inv={matching.inv})")
    
    # Show examples where inversion rule fails
    inv_fails = {k: v for k, v in all_valid_by_inv.items() 
                 if len([t for t in v if t[0] == 'inv_match']) != 1}
    
    if inv_fails:
        print()
        print(f"✗ Inversion rule FAILS for {len(inv_fails)} pairs")
        print(f"  Examples (first 3):")
        for idx, ((rc1, rc2), terms) in enumerate(list(inv_fails.items())[:3]):
            matches = [t[1] for t in terms if t[0] == 'inv_match']
            print(f"    Pair ({rc1.perm.trimcode}, {rc2.perm.trimcode})")
            print(f"      Inversions: {rc1.perm.inv} + {rc2.perm.inv} = {rc1.perm.inv + rc2.perm.inv}")
            print(f"      Matches {len(matches)} sources: {[m.trimcode for m in matches]}")
    
    print()
    
    # Test Proposal D: Maximal element
    print("="*60)
    print("TESTING PROPOSAL D: HIGHEST WEIGHT CONDITION")
    print("="*60)
    print()
    
    # For each (perm1, perm2) pair across ALL coproducts,
    # check if choosing the highest weight representative gives uniqueness
    
    all_by_perm_pair = {}  # Maps (perm1, perm2) -> list of (rc1, rc2, source_perm)
    
    for perm in perms:
        principal = RCGraph.principal_rc(perm, n-1)
        coprod = rc_ring.coproduct_on_basis(principal)
        
        for (rc1, rc2), coeff in coprod.items():
            perm_pair = (rc1.perm, rc2.perm)
            if perm_pair not in all_by_perm_pair:
                all_by_perm_pair[perm_pair] = []
            all_by_perm_pair[perm_pair].append((rc1, rc2, perm))
    
    # For each perm pair, check if highest weight choice is unique
    hw_unique_count = 0
    hw_ambiguous = []
    
    for perm_pair, rc_list in all_by_perm_pair.items():
        # Get highest weight for each RC pair
        hw_data = []
        for rc1, rc2, src in rc_list:
            hw1, _ = rc1.to_highest_weight()
            hw2, _ = rc2.to_highest_weight()
            hw_data.append((tuple(hw1.crystal_weight), tuple(hw2.crystal_weight), rc1, rc2, src))
        
        # Find maximum
        max_hw = max(hw_data, key=lambda x: (x[0], x[1]))
        max_hw_count = sum(1 for hwd in hw_data if (hwd[0], hwd[1]) == (max_hw[0], max_hw[1]))
        
        if max_hw_count == 1:
            hw_unique_count += 1
        else:
            hw_ambiguous.append((perm_pair, max_hw_count, [hwd[4] for hwd in hw_data if (hwd[0], hwd[1]) == (max_hw[0], max_hw[1])]))
    
    print(f"Permutation pairs with UNIQUE highest weight: {hw_unique_count} / {len(all_by_perm_pair)}")
    
    if hw_ambiguous:
        print(f"✗ Highest weight ambiguous for {len(hw_ambiguous)} pairs")
        print(f"  Examples (first 3):")
        for idx, (pp, count, sources) in enumerate(hw_ambiguous[:3]):
            print(f"    Pair ({pp[0].trimcode}, {pp[1].trimcode})")
            print(f"      {count} RC pairs with same highest weight")
            print(f"      From sources: {[s.trimcode for s in sources]}")
    
    print()


if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    
    test_local_canonicality_conditions(n)
    propose_local_canonical_condition(n)