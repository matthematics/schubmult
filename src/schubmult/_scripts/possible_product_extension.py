"""
Extend the RC-graph product using the twisted Leibniz rule.

Strategy: Use the Leibniz formula recursively to define products
when the Leibniz rule forces a unique answer.
"""

def analyze_leibniz_extension():
    """
    The twisted Leibniz rule is:
        ∂_i(a % b) = ∂_i(a) % b + s_i(a) % ∂_i(b)
    
    Rearranging:
        a % b = ???  (if we know ∂_i(a % b), ∂_i(a), s_i(a), ∂_i(b))
    
    This gives several extension strategies:
    
    STRATEGY 1: Extend when ∂_i(left) = 0
    ----------------------------------------
    If ∂_i(a) = 0, then:
        ∂_i(a % b) = s_i(a) % ∂_i(b)
    
    This DEFINES a % ∂_i(b) in terms of known products!
    
    Algorithm:
    - If we know a % dom (where dom is dominant)
    - And ∂_i(a) = 0
    - Then we can define: a % ∂_i(dom) := ∂_i(a % dom)
    
    This extends the product to NON-DOMINANT right factors!
    
    
    STRATEGY 2: Extend when s_i(a) is simpler
    ------------------------------------------
    If we know:
    - ∂_i(a % b) 
    - ∂_i(a) % b
    - ∂_i(b)
    
    Then: s_i(a) % ∂_i(b) = ∂_i(a % b) - ∂_i(a) % b
    
    This defines s_i(a) % ∂_i(b) if s_i(a) is "new".
    
    
    STRATEGY 3: Lifting via highest weight
    ---------------------------------------
    Every RC-graph rc can be written as:
        rc = hw.reverse_raise_seq(seq)
    
    where hw is highest weight.
    
    If we know hw1 % hw2, we can compute:
        rc1 % rc2 by working in crystal components
    """
    print("\n" + "="*80)
    print("STRATEGIES FOR EXTENDING % PRODUCT VIA LEIBNIZ RULE")
    print("="*80)
    
    print("\nCurrent status:")
    print("  ✓ dom % any  (dominant left factor)")
    print("  ✓ any % dom  (dominant right factor)")
    print("  ? any % any  (general case)")
    print()
    
    print("Extension strategies:")
    print()
    print("1. USE ∂_i(left) = 0:")
    print("   If ∂_i(a) = 0 and we know a % dom:")
    print("     Define a % ∂_i(dom) := ∂_i(a % dom)")
    print()
    print("2. USE LEIBNIZ RECURSIVELY:")
    print("   If ∂_i(a % b), ∂_i(a) % b are known:")
    print("     Define s_i(a) % ∂_i(b) := ∂_i(a % b) - ∂_i(a) % b")
    print()
    print("3. USE CRYSTAL STRUCTURE:")
    print("   Express both factors in terms of highest weights + lowering")
    print("   Compute product of highest weights, then lower")
    print()


def check_leibniz_for_dominant_right(n, verbose=True):
    """
    Verify that the twisted Leibniz rule holds for % product
    when right factor is dominant.
    
    Tests: ∂_i(a % dom) = ∂_i(a) % dom + s_i(a) % ∂_i(dom)
    
    This should hold by construction if rc_single_product is correct.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"VERIFYING LEIBNIZ RULE FOR % PRODUCT (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant()]
    
    rc_ring = RCGraphRing()
    
    passed = 0
    failed = 0
    errors = 0
    
    for a_perm in perms[:5]:  # Sample
        for dom_perm in dom_perms:
            for i in range(1, n):
                try:
                    a_rc = RCGraph.principal_rc(a_perm, n-1)
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    
                    # LHS: ∂_i(a % dom)
                    prod = rc_ring(a_rc) % rc_ring(dom_rc)
                    lhs = prod.divdiff(i)
                    
                    # RHS term 1: ∂_i(a) % dom
                    diff_a = rc_ring(a_rc).divdiff(i)
                    term1 = diff_a % rc_ring(dom_rc)
                    
                    # RHS term 2: s_i(a) % ∂_i(dom)
                    s_i = Permutation.simple_transposition(n, i)
                    
                    # Apply s_i to a_rc (need to implement this!)
                    # For now, approximate via divided difference
                    if i in a_perm.descents():
                        s_i_a = diff_a  # s_i acts as ∂_i on descents
                    else:
                        # s_i multiplies on left
                        s_i_a_perm = s_i * a_perm
                        s_i_a = rc_ring(RCGraph.principal_rc(s_i_a_perm, n-1))
                    
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    term2 = s_i_a % diff_dom
                    
                    rhs = term1 + term2
                    
                    # Compare
                    diff = lhs - rhs
                    
                    if all(abs(v) < 1e-10 for v in diff.values()):
                        passed += 1
                        if verbose:
                            print(f"  ✓ Leibniz holds: a={a_perm.trimcode}, dom={dom_perm.trimcode}, i={i}")
                    else:
                        failed += 1
                        print(f"\n  ✗ Leibniz FAILS: a={a_perm.trimcode}, dom={dom_perm.trimcode}, i={i}")
                        print(f"    LHS = {lhs}")
                        print(f"    RHS = {rhs}")
                        print(f"    Diff = {diff}")
                
                except Exception as e:
                    errors += 1
                    if verbose:
                        print(f"  ⚠ Error: a={a_perm.trimcode}, dom={dom_perm.trimcode}, i={i}")
                        print(f"    {e}")
    
    print(f"\n=== Summary ===")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"Errors: {errors}")
    
    if failed == 0 and errors == 0:
        print("\n✓ Leibniz rule VERIFIED for dominant right factor!")
    
    return passed, failed, errors


def extend_product_via_zero_divdiff(rc_ring, a_rc, dom_rc, i):
    """
    Strategy 1: Extend product when ∂_i(a) = 0.
    
    If ∂_i(a) = 0 and we know a % dom, then:
        a % ∂_i(dom) := ∂_i(a % dom)
    
    This extends the product to non-dominant right factors.
    """
    # Check if ∂_i(a) = 0
    diff_a = rc_ring(a_rc).divdiff(i)
    
    if len(diff_a) == 0:  # ∂_i(a) = 0
        # Compute a % dom (known, since dom is dominant)
        prod = rc_ring(a_rc) % rc_ring(dom_rc)
        
        # Define a % ∂_i(dom) := ∂_i(a % dom)
        diff_dom = rc_ring(dom_rc).divdiff(i)
        
        # For each rc in diff_dom, we've now defined a % rc
        new_products = {}
        prod_diff = prod.divdiff(i)
        
        for rc_b, coeff_b in diff_dom.items():
            # a % rc_b should equal (1/coeff_b) * prod_diff[corresponding terms]
            # This is approximate - need to match crystal components
            new_products[rc_b] = prod_diff  # Simplified
        
        return new_products
    
    return None


def find_zero_divdiff_pairs(n):
    """
    Find all pairs (a, i) where ∂_i(a) = 0.
    
    These are the cases where Strategy 1 applies.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"FINDING ZERO DIVIDED DIFFERENCE PAIRS (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    rc_ring = RCGraphRing()
    
    zero_pairs = []
    
    for perm in perms:
        rc = RCGraph.principal_rc(perm, n-1)
        
        for i in range(1, n):
            diff = rc_ring(rc).divdiff(i)
            
            if len(diff) == 0:
                zero_pairs.append((perm, i))
    
    print(f"\nFound {len(zero_pairs)} pairs (a, i) with ∂_i(a) = 0:")
    print()
    
    # Group by permutation
    by_perm = {}
    for perm, i in zero_pairs:
        if perm not in by_perm:
            by_perm[perm] = []
        by_perm[perm].append(i)
    
    for perm in sorted(by_perm.keys(), key=lambda p: p.inv):
        indices = by_perm[perm]
        print(f"  {perm.trimcode}: ∂_i = 0 for i ∈ {{{', '.join(map(str, indices))}}}")
    
    print("\nPattern analysis:")
    print("  These are permutations where certain positions are 'stable'")
    print("  under the divided difference operator.")
    print("  For these pairs, we can extend the product to non-dominant right factors!")
    
    return zero_pairs


def propose_systematic_extension():
    """
    Propose a systematic algorithm to extend % product to all pairs.
    """
    print("\n" + "="*80)
    print("SYSTEMATIC EXTENSION ALGORITHM")
    print("="*80)
    
    print("\nAlgorithm: Extend % product to all RC-graph pairs")
    print()
    print("INPUT: Two RC-graphs a_rc, b_rc")
    print("OUTPUT: a_rc % b_rc (as RCGraphRingElement)")
    print()
    print("STEP 1: Check if b_rc is dominant")
    print("  If YES: Use existing rc_single_product")
    print("  If NO: Proceed to Step 2")
    print()
    print("STEP 2: Express b_rc in terms of dominant + divided differences")
    print("  b_rc = ∂_{i_k}...∂_{i_1}(dom_rc)")
    print("  where dom_rc = minimal_dominant_above(b_rc.perm)")
    print()
    print("STEP 3: Apply Leibniz rule recursively")
    print("  For each j from 1 to k:")
    print("    Use: ∂_{i_j}(a % b') = ∂_{i_j}(a) % b' + s_{i_j}(a) % ∂_{i_j}(b')")
    print()
    print("STEP 4: Handle base cases")
    print("  - If ∂_i(a) = 0: Use Strategy 1 (extend directly)")
    print("  - If s_i(a) is known: Use Strategy 2 (recursive definition)")
    print("  - Otherwise: Use crystal decomposition")
    print()
    print("STEP 5: Consistency checks")
    print("  - Verify result is independent of choice of divided difference sequence")
    print("  - Check associativity: (a % b) % c = a % (b % c)")
    print("  - Verify polynomial equality: (a % b).polyvalue() = Sx(a.perm) * Sx(b.perm)")
    print()


def implement_extended_product():
    """
    Implement extended % product using Leibniz rule.
    
    Add this to RCGraphRing class.
    """
    code = '''
def rc_single_product_extended(self, u_rc, v_rc):
    """
    Extended % product using Leibniz rule.
    
    Handles:
    - v_rc dominant: Use existing implementation
    - v_rc non-dominant: Use Leibniz rule recursively
    """
    # Base case: v_rc is dominant
    if v_rc.perm.is_dominant():
        return self.rc_single_product(u_rc, v_rc)
    
    # Find divided difference representation of v_rc
    dom_v = v_rc.perm.minimal_dominant_above()
    dom_v_rc = RCGraph.principal_rc(dom_v, len(v_rc))
    
    # Find sequence of divided differences: v_rc = ∂_{i_k}...∂_{i_1}(dom_v_rc)
    diff_seq = find_divdiff_sequence(v_rc, dom_v_rc)
    
    # Apply Leibniz rule recursively
    # Start with u_rc % dom_v_rc (known)
    result = self.rc_single_product(u_rc, dom_v_rc)
    
    # Apply each divided difference in reverse
    for i in reversed(diff_seq):
        # ∂_i(u_rc % current) = ∂_i(u_rc) % current + s_i(u_rc) % ∂_i(current)
        
        # Check if ∂_i(u_rc) = 0 (Strategy 1)
        diff_u = self(u_rc).divdiff(i)
        
        if len(diff_u) == 0:
            # ∂_i(u_rc) = 0, so just apply ∂_i to result
            result = result.divdiff(i)
        else:
            # General case: need to implement s_i action
            # This is the tricky part!
            raise NotImplementedError("Need to implement s_i action on RC-graphs")
    
    return result


def find_divdiff_sequence(rc, dom_rc):
    """
    Find sequence of divided differences connecting rc to dom_rc.
    
    Returns list [i_1, ..., i_k] such that:
        rc = ∂_{i_k}...∂_{i_1}(dom_rc)
    """
    # Algorithm: Walk from dom_rc down to rc via lowering operators
    current = dom_rc
    seq = []
    
    while current.perm != rc.perm:
        # Find which lowering operator to apply
        for i in range(1, len(rc) + 1):
            lowered = current.lowering_operator(i)
            if lowered is not None and lowered.perm.bruhat_leq(rc.perm):
                seq.append(i)
                current = lowered
                break
        else:
            raise ValueError("Cannot find divided difference sequence")
    
    return seq
'''
    
    print("\n" + "="*80)
    print("IMPLEMENTATION: EXTENDED % PRODUCT")
    print("="*80)
    print()
    print("Add to RCGraphRing class:")
    print()
    print(code)
    print()
    print("Key challenge: Implementing s_i action on RC-graphs")
    print("  This requires understanding Demazure operator structure")
    print()


"""
Refined understanding of divided differences and crystal structure.

Key correction: divdiff_desc(i) ≠ 0 requires TWO conditions:
1. i is a descent of rc.perm
2. The simple root α_i appears in row i of the RC-graph

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
    print()


def when_is_divdiff_zero(n):
    """
    Analyze when divdiff(rc, i) = 0.
    
    This happens when BOTH:
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
    
    for perm in perms:
        rc = RCGraph.principal_rc(perm, n-1)
        
        for i in range(1, n):
            diff = rc_ring(rc).divdiff(i)
            
            if len(diff) == 0:
                # Classify why it's zero
                if i not in perm.descents():
                    zero_cases['no_descent'].append((perm, i))
                else:
                    # i is a descent, check if α_i is in row i
                    # This requires inspecting the RC-graph structure
                    desc_result = rc.divdiff_desc(i)
                    
                    if desc_result is None:
                        # α_i not in row i, or doesn't yield hw
                        zero_cases['no_root'].append((perm, i))
                    else:
                        # Shouldn't happen if divdiff = 0 but divdiff_desc ≠ 0
                        print(f"  ⚠ Unexpected: {perm.trimcode}, i={i}")
                        print(f"    divdiff = 0 but divdiff_desc ≠ 0")
    
    print(f"\nCases where divdiff(rc, i) = 0:")
    print()
    print(f"1. No descent at i: {len(zero_cases['no_descent'])} cases")
    print(f"2. Descent but no α_i in row i: {len(zero_cases['no_root'])} cases")
    print()
    
    # Show examples
    print("Examples of 'no descent' cases:")
    for perm, i in zero_cases['no_descent'][:5]:
        print(f"  {perm.trimcode}: i={i} (not in {perm.descents()})")
    
    print("\nExamples of 'descent but no root' cases:")
    for perm, i in zero_cases['no_root'][:5]:
        print(f"  {perm.trimcode}: i={i} (in {perm.descents()} but α_i not in row i)")
    
    return zero_cases


def find_useful_zero_divdiff_pairs(n):
    """
    Find pairs (a, i) where divdiff(a, i) = 0 AND this is useful
    for extending products via Leibniz rule.
    
    These are cases where we can use:
        ∂_i(a % dom) = s_i(a) % ∂_i(dom)
    
    to extend the product to non-dominant right factors.
    """
    from schubmult import Permutation, RCGraph, RCGraphRing
    
    print("\n" + "="*80)
    print(f"USEFUL ZERO DIVDIFF PAIRS FOR EXTENSION (n={n})")
    print("="*80)
    
    perms = Permutation.all_permutations(n)
    dom_perms = [p for p in perms if p.is_dominant()]
    rc_ring = RCGraphRing()
    
    useful_pairs = []
    
    for a_perm in perms:
        a_rc = RCGraph.principal_rc(a_perm, n-1)
        
        for i in range(1, n):
            diff_a = rc_ring(a_rc).divdiff(i)
            
            if len(diff_a) == 0:
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
                    useful_pairs.append((a_perm, i))
    
    print(f"\nFound {len(useful_pairs)} useful pairs where:")
    print("  - divdiff(a, i) = 0")
    print("  - ∃ dominant dom with divdiff(dom, i) ≠ 0")
    print()
    
    # Group by permutation
    by_perm = {}
    for perm, i in useful_pairs:
        if perm not in by_perm:
            by_perm[perm] = []
        by_perm[perm].append(i)
    
    print("Useful pairs by permutation:")
    for perm in sorted(by_perm.keys(), key=lambda p: p.inv)[:10]:
        indices = by_perm[perm]
        print(f"  {perm.trimcode}: i ∈ {{{', '.join(map(str, indices))}}}")
    
    print("\nFor these pairs, we can extend:")
    print("  a % ∂_i(dom)  is determined by  ∂_i(a % dom)")
    print("  via the Leibniz formula:")
    print("    ∂_i(a % dom) = s_i(a) % ∂_i(dom)")
    print()
    
    return useful_pairs


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
    
    for perm in perms[:5]:
        rc = RCGraph.principal_rc(perm, n-1)
        
        for i in range(1, n):
            s_i_perm = Permutation.simple_transposition(n, i)
            
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
    dom_perms = [p for p in perms if p.is_dominant()]
    rc_ring = RCGraphRing()
    
    # Find a test case
    test_cases = []
    
    for a_perm in perms[:10]:
        a_rc = RCGraph.principal_rc(a_perm, n-1)
        
        for i in range(1, n):
            diff_a = rc_ring(a_rc).divdiff(i)
            
            if len(diff_a) == 0:
                # Found a useful case
                for dom_perm in dom_perms[:2]:
                    dom_rc = RCGraph.principal_rc(dom_perm, n-1)
                    diff_dom = rc_ring(dom_rc).divdiff(i)
                    
                    if len(diff_dom) > 0:
                        test_cases.append((a_rc, dom_rc, i))
                        break
                break
        
        if len(test_cases) >= 3:
            break
    
    print(f"\nTesting {len(test_cases)} cases:")
    print()
    
    for a_rc, dom_rc, i in test_cases:
        print(f"Test: a={a_rc.perm.trimcode}, dom={dom_rc.perm.trimcode}, i={i}")
        
        # LHS: ∂_i(a % dom)
        prod_a_dom = rc_ring(a_rc) % rc_ring(dom_rc)
        lhs = prod_a_dom.divdiff(i)
        
        print(f"  LHS = ∂_{i}(a % dom) has {len(lhs)} terms")
        
        # RHS: s_i(a) % ∂_i(dom)
        
        # Compute s_i(a)
        s_i_perm = Permutation.simple_transposition(n, i)
        
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
        
        # The equation LHS = RHS gives us the information
        # to define each product individually
        
        # Match by crystal weights
        lhs_weights = {}
        for rc, coeff in lhs.items():
            weight = tuple(rc.crystal_weight)
            if weight not in lhs_weights:
                lhs_weights[weight] = []
            lhs_weights[weight].append((rc, coeff))
        
        print(f"  LHS has {len(lhs_weights)} distinct weights")
        
        # For each product s_i(a_j) % dom_k, predict the weight
        for s_i_a_rc, s_i_a_coeff in s_i_a.items():
            for dom_k, dom_k_coeff in diff_dom.items():
                # Expected weight
                w_a = s_i_a_rc.crystal_weight
                w_d = dom_k.crystal_weight
                expected_weight = tuple(w_a[j] + w_d[j] for j in range(len(w_a)))
                
                print(f"    Product {s_i_a_rc.perm.trimcode} % {dom_k.perm.trimcode}")
                print(f"      Expected weight: {expected_weight}")
                
                if expected_weight in lhs_weights:
                    print(f"      Found in LHS! ({len(lhs_weights[expected_weight])} terms)")
                else:
                    print(f"      NOT in LHS - unexpected!")
        
        print()


if __name__ == "__main__":
    from schubmult import Permutation, RCGraph, RCGraphRing, Sx

    import sys
    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    
    clarify_divdiff_structure()
    
    zero_cases = when_is_divdiff_zero(n)
    
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
    print("For Leibniz extension when divdiff(a, i) = 0:")
    print("  ∂_i(a % dom) = s_i(a) % ∂_i(dom)")
    print()
    print("Where s_i(a) is:")
    print("  - ∂_i(a) if i is a descent")
    print("  - s_i · a otherwise")
    print()
    print(f"Found {len(useful_pairs)} useful pairs for extension!")

    
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    
    analyze_leibniz_extension()
    
    # Check that Leibniz holds for current implementation
    check_leibniz_for_dominant_right(n, verbose=False)
    
    # Find cases where Strategy 1 applies
    zero_pairs = find_zero_divdiff_pairs(n)
    
    # Propose systematic extension
    propose_systematic_extension()
    
    # Show implementation
    implement_extended_product()
    
    print("\n" + "="*80)
    print("SUMMARY: EXTENDING % PRODUCT")
    print("="*80)
    print()
    print("Current status:")
    print("  ✓ any % dom  (dominant right factor)")
    print()
    print("Extension strategies:")
    print("  1. ∂_i(left) = 0  →  Extend to ∂_i(dom) right factors")
    print("  2. Leibniz recursion  →  General extension")
    print("  3. Crystal structure  →  Highest weight decomposition")
    print()
    print("Next steps:")
    print("  1. Verify Leibniz rule holds for current implementation")
    print(f"  2. Identify {len(zero_pairs)} zero-divdiff pairs for Strategy 1")
    print("  3. Implement s_i action on RC-graphs")
    print("  4. Extend product recursively via Leibniz formula")
    print()
    print("Key insight: The twisted Leibniz rule FORCES the extension!")
    print("There is (likely) a UNIQUE way to extend % to all pairs.")



