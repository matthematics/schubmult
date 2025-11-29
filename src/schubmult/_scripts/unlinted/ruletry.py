# Add to more_claude.py

def analyze_coproduct_cancellation_patterns(n):
    """
    Analyze the cancellation patterns in RCGraphRing.coproduct().
    
    Goal: Find a canonical way to identify which RC-graphs cancel,
    beyond just matching permutations.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"ANALYZING COPRODUCT CANCELLATION PATTERNS (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    print("\nGoal: Find canonical parametrization to identify which RC-graphs")
    print("      cancel in the coproduct, beyond just matching permutations.")
    print()
    
    # Get some RC-graphs to test
    test_rcs = []
    for perm in perms[:10]:  # Test first 10 perms
        principal = RCGraph.principal_rc(perm, n-1)
        test_rcs.append(principal)
        
        # Add some non-principal if they exist
        all_rcs = RCGraph.all_rc_graphs(perm, n-1)
        for rc in all_rcs[:2]:  # Up to 2 per perm
            if rc != principal:
                test_rcs.append(rc)
    
    print(f"Testing {len(test_rcs)} RC-graphs")
    print()
    
    # For each RC-graph, compute its coproduct
    cancellation_data = []
    
    for idx, rc in enumerate(test_rcs):
        print("="*60)
        print(f"RC-GRAPH {idx+1}/{len(test_rcs)}: perm {rc.perm.trimcode}")
        print(f"  Weight: {rc.crystal_weight}")
        
        principal = RCGraph.principal_rc(rc.perm, n-1)
        is_principal = (rc == principal)
        print(f"  Is principal? {is_principal}")
        
        pretty_print(rc)
        print()
        
        # Compute coproduct
        coprod = rc_ring(rc).coproduct()
        
        print(f"Coproduct has {len(coprod)} terms (before cancellation)")
        
        # Group by permutation pairs
        by_perm_pair = {}
        for (rc1, rc2), coeff in coprod.items():
            perm_pair = (rc1.perm.trimcode, rc2.perm.trimcode)
            if perm_pair not in by_perm_pair:
                by_perm_pair[perm_pair] = []
            by_perm_pair[perm_pair].append(((rc1, rc2), coeff))
        
        print(f"Grouped into {len(by_perm_pair)} distinct permutation pairs")
        print()
        
        # Analyze each permutation pair that has potential for cancellation
        for perm_pair, terms in by_perm_pair.items():
            if len(terms) <= 1:
                continue  # No cancellation possible
            
            print(f"Permutation pair {perm_pair}:")
            print(f"  {len(terms)} terms with this permutation pair")
            print()
            
            # For each term, compute various canonical identifiers
            for (rc1, rc2), coeff in terms:
                print(f"  Term: coeff={coeff}")
                print(f"    Left RC:")
                print(f"      Weight: {rc1.crystal_weight}")
                
                # Canonical ID 1: Highest weight + lowering sequence
                hw1, raise_seq1 = rc1.to_highest_weight()
                lowering_seq1 = tuple(reversed(raise_seq1))
                hw_weight1 = tuple(hw1.crystal_weight)
                
                print(f"      HW weight: {hw_weight1}")
                print(f"      Lowering seq: {lowering_seq1}")
                
                # Canonical ID 2: Reading word of roots
                # (Read roots row by row, left to right)
                reading_word1 = []
                for i, row in enumerate(rc1):
                    for root in sorted(row):
                        reading_word1.append((i, root))
                reading_word1 = tuple(reading_word1)
                print(f"      Reading word: {reading_word1}")
                
                # Canonical ID 3: Params (if available)
                try:
                    params1 = rc1.params
                    print(f"      Params: {params1}")
                except:
                    params1 = None
                    print(f"      Params: N/A")
                
                print(f"    Right RC:")
                print(f"      Weight: {rc2.crystal_weight}")
                
                hw2, raise_seq2 = rc2.to_highest_weight()
                lowering_seq2 = tuple(reversed(raise_seq2))
                hw_weight2 = tuple(hw2.crystal_weight)
                
                print(f"      HW weight: {hw_weight2}")
                print(f"      Lowering seq: {lowering_seq2}")
                
                reading_word2 = []
                for i, row in enumerate(rc2):
                    for root in sorted(row):
                        reading_word2.append((i, root))
                reading_word2 = tuple(reading_word2)
                print(f"      Reading word: {reading_word2}")
                
                try:
                    params2 = rc2.params
                    print(f"      Params: {params2}")
                except:
                    params2 = None
                    print(f"      Params: N/A")
                
                # Store all identifiers
                cancellation_data.append({
                    'original_rc': rc.perm.trimcode,
                    'perm_pair': perm_pair,
                    'coeff': coeff,
                    'left': {
                        'weight': rc1.crystal_weight,
                        'hw_weight': hw_weight1,
                        'lowering_seq': lowering_seq1,
                        'reading_word': reading_word1,
                        'params': params1,
                    },
                    'right': {
                        'weight': rc2.crystal_weight,
                        'hw_weight': hw_weight2,
                        'lowering_seq': lowering_seq2,
                        'reading_word': reading_word2,
                        'params': params2,
                    }
                })
                print()
            
            # Check if any identifiers distinguish the terms
            print(f"  Analysis of distinguishing features:")
            
            # Group by different canonical IDs
            by_hw_and_lowering = {}
            by_reading_word = {}
            by_params = {}
            
            for (rc1, rc2), coeff in terms:
                hw1, raise_seq1 = rc1.to_highest_weight()
                hw2, raise_seq2 = rc2.to_highest_weight()
                
                id_hw_low = (
                    tuple(hw1.crystal_weight), tuple(reversed(raise_seq1)),
                    tuple(hw2.crystal_weight), tuple(reversed(raise_seq2))
                )
                
                rw1 = tuple((i, r) for i, row in enumerate(rc1) for r in sorted(row))
                rw2 = tuple((i, r) for i, row in enumerate(rc2) for r in sorted(row))
                id_reading = (rw1, rw2)
                
                try:
                    id_params = (rc1.params, rc2.params)
                except:
                    id_params = None
                
                if id_hw_low not in by_hw_and_lowering:
                    by_hw_and_lowering[id_hw_low] = []
                by_hw_and_lowering[id_hw_low].append(coeff)
                
                if id_reading not in by_reading_word:
                    by_reading_word[id_reading] = []
                by_reading_word[id_reading].append(coeff)
                
                if id_params is not None:
                    if id_params not in by_params:
                        by_params[id_params] = []
                    by_params[id_params].append(coeff)
            
            print(f"    By (HW weight + lowering seq): {len(by_hw_and_lowering)} distinct values")
            for id_val, coeffs in by_hw_and_lowering.items():
                print(f"      {id_val}: {coeffs} (sum={sum(coeffs)})")
            
            print(f"    By reading word: {len(by_reading_word)} distinct values")
            for id_val, coeffs in by_reading_word.items():
                print(f"      {id_val}: {coeffs} (sum={sum(coeffs)})")
            
            if by_params:
                print(f"    By params: {len(by_params)} distinct values")
                for id_val, coeffs in by_params.items():
                    print(f"      {id_val}: {coeffs} (sum={sum(coeffs)})")
            
            print()
    
    return cancellation_data


def test_canonical_cancellation_methods(n):
    """
    Test different canonical cancellation methods to see which gives
    the most natural positive formula.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"TESTING CANONICAL CANCELLATION METHODS (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    print("\nWe will test different canonical identifiers:")
    print("  Method 1: (HW weight, lowering sequence)")
    print("  Method 2: Reading word of roots")
    print("  Method 3: Params (if available)")
    print()
    
    # Test on several RC-graphs
    test_perms = perms[:5]
    
    for perm in test_perms:
        print("="*60)
        print(f"Testing perm {perm.trimcode}")
        
        principal = RCGraph.principal_rc(perm, n-1)
        all_rcs = RCGraph.all_rc_graphs(perm, n-1)
        
        print(f"  Total RC-graphs: {len(all_rcs)}")
        print()
        
        for rc in all_rcs[:3]:  # Test first 3
            is_principal = (rc == principal)
            marker = " (principal)" if is_principal else ""
            
            print(f"RC-graph{marker}:")
            print(f"  Weight: {rc.crystal_weight}")
            pretty_print(rc)
            print()
            
            # Get coproduct
            coprod = rc_ring(rc).coproduct()
            
            print(f"Original coproduct: {len(coprod)} terms")
            
            # Method 1: Cancel by (HW weight, lowering seq)
            method1_cancelled = cancel_by_hw_lowering(coprod)
            print(f"Method 1 (HW + lowering): {len(method1_cancelled)} terms after cancellation")
            
            # Check if all coefficients are positive
            all_positive_1 = all(coeff > 0 for coeff in method1_cancelled.values())
            print(f"  All coefficients positive? {all_positive_1}")
            if not all_positive_1:
                print(f"  Negative coefficients found!")
                for (rc1, rc2), coeff in method1_cancelled.items():
                    if coeff < 0:
                        print(f"    ({rc1.perm.trimcode}, {rc2.perm.trimcode}): {coeff}")
            
            # Method 2: Cancel by reading word
            method2_cancelled = cancel_by_reading_word(coprod)
            print(f"Method 2 (reading word): {len(method2_cancelled)} terms after cancellation")
            
            all_positive_2 = all(coeff > 0 for coeff in method2_cancelled.values())
            print(f"  All coefficients positive? {all_positive_2}")
            if not all_positive_2:
                print(f"  Negative coefficients found!")
            
            # Method 3: Cancel by params (if available)
            try:
                method3_cancelled = cancel_by_params(coprod)
                print(f"Method 3 (params): {len(method3_cancelled)} terms after cancellation")
                
                all_positive_3 = all(coeff > 0 for coeff in method3_cancelled.values())
                print(f"  All coefficients positive? {all_positive_3}")
                if not all_positive_3:
                    print(f"  Negative coefficients found!")
            except Exception as e:
                print(f"Method 3 (params): Failed - {e}")
            
            print()


def cancel_by_hw_lowering(coprod_dict):
    """
    Cancel terms in coproduct by (HW weight, lowering sequence) ID.
    """
    by_canonical = {}
    
    for (rc1, rc2), coeff in coprod_dict.items():
        hw1, raise_seq1 = rc1.to_highest_weight()
        hw2, raise_seq2 = rc2.to_highest_weight()
        
        canonical_id = (
            rc1.perm.trimcode, rc2.perm.trimcode,
            tuple(hw1.crystal_weight), tuple(reversed(raise_seq1)),
            tuple(hw2.crystal_weight), tuple(reversed(raise_seq2))
        )
        
        if canonical_id not in by_canonical:
            by_canonical[canonical_id] = []
        by_canonical[canonical_id].append(((rc1, rc2), coeff))
    
    # Sum coefficients for each canonical ID
    result = {}
    for canonical_id, terms in by_canonical.items():
        total_coeff = sum(c for _, c in terms)
        if total_coeff != 0:
            # Take first representative
            result[terms[0][0]] = total_coeff
    
    return result


def cancel_by_reading_word(coprod_dict):
    """
    Cancel terms in coproduct by reading word ID.
    """
    by_canonical = {}
    
    for (rc1, rc2), coeff in coprod_dict.items():
        rw1 = tuple((i, r) for i, row in enumerate(rc1) for r in sorted(row))
        rw2 = tuple((i, r) for i, row in enumerate(rc2) for r in sorted(row))
        
        canonical_id = (rc1.perm.trimcode, rc2.perm.trimcode, rw1, rw2)
        
        if canonical_id not in by_canonical:
            by_canonical[canonical_id] = []
        by_canonical[canonical_id].append(((rc1, rc2), coeff))
    
    result = {}
    for canonical_id, terms in by_canonical.items():
        total_coeff = sum(c for _, c in terms)
        if total_coeff != 0:
            result[terms[0][0]] = total_coeff
    
    return result


def cancel_by_params(coprod_dict):
    """
    Cancel terms in coproduct by params.
    """
    by_canonical = {}
    
    for (rc1, rc2), coeff in coprod_dict.items():
        canonical_id = (rc1.perm.trimcode, rc2.perm.trimcode, rc1.params, rc2.params)
        
        if canonical_id not in by_canonical:
            by_canonical[canonical_id] = []
        by_canonical[canonical_id].append(((rc1, rc2), coeff))
    
    result = {}
    for canonical_id, terms in by_canonical.items():
        total_coeff = sum(c for _, c in terms)
        if total_coeff != 0:
            result[terms[0][0]] = total_coeff
    
    return result


if __name__ == "__main__":
    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    
    # ... existing functions ...
    
    # NEW: Analyze coproduct cancellation patterns
    cancellation_data = analyze_coproduct_cancellation_patterns(n)
    
    # NEW: Test different canonical cancellation methods
    test_canonical_cancellation_methods(n)
    
    # ... rest of existing code ...