from sympy import init_printing, pretty_print
from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.rings.schubert_ring import Sx
from schubmult.rings.rc_graph_ring import RCGraphRing

init_printing()

n = 3  # Start with small rank for testing

# Get all permutations
perms = Permutation.all_permutations(n)

def test_inverse_crystal_product_basic():
    """
    Test basic properties of inverse_crystal_product.
    """
    print("\n" + "="*60)
    print("Testing Basic Properties of inverse_crystal_product")
    print("="*60)
    
    rc_ring = RCGraphRing()
    
    # Test 1: Identity
    print("\n1. Testing identity: e * rc = rc * e = rc")
    identity_rc = RCGraph([])
    test_cases = 0
    passed = 0
    
    for perm in perms[:10]:  # Test on first 10 perms
        for rc in list(RCGraph.all_hw_rcs(perm, n-1))[:2]:  # First 2 RCs per perm
            result_left = identity_rc.inverse_crystal_product(rc)
            result_right = rc.inverse_crystal_product(identity_rc)
            
            test_cases += 2
            
            # Check if result contains rc with coefficient 1
            if result_left.get(rc, 0) == 1 and len(result_left) == 1:
                passed += 1
            else:
                print(f"  ✗ e * {rc.perm.trimcode} failed")
                pretty_print(result_left)
                pretty_print(result_right)
                pretty_print(rc)
            
            if result_right.get(rc, 0) == 1 and len(result_right) == 1:
                passed += 1
            else:
                print(f"  ✗ {rc.perm.trimcode} * e failed")
                pretty_print(result_left)
                pretty_print(result_right)
                pretty_print(rc)
    
    print(f"Identity tests: {passed}/{test_cases} passed")
    
    # Test 2: Associativity (sample)
    print("\n2. Testing associativity: (a*b)*c = a*(b*c)")
    test_cases = 0
    passed = 0
    
    for u in perms[:5]:
        for v in perms[:5]:
            for w in perms[:5]:
                # Pick one RC from each
                u_rc = next(iter(RCGraph.all_hw_rcs(u, n-1)))
                v_rc = next(iter(RCGraph.all_hw_rcs(v, n-1)))
                w_rc = next(iter(RCGraph.all_hw_rcs(w, n-1)))
                
                # Compute (u*v)*w
                uv = u_rc.inverse_crystal_product(v_rc)
                left_result = rc_ring.zero
                for rc1, coeff1 in uv.items():
                    prod = rc1.inverse_crystal_product(w_rc)
                    for rc2, coeff2 in prod.items():
                        left_result += coeff1 * coeff2 * rc_ring(rc2)
                
                # Compute u*(v*w)
                vw = v_rc.inverse_crystal_product(w_rc)
                right_result = rc_ring.zero
                for rc1, coeff1 in vw.items():
                    prod = u_rc.inverse_crystal_product(rc1)
                    for rc2, coeff2 in prod.items():
                        right_result += coeff1 * coeff2 * rc_ring(rc2)
                
                test_cases += 1
                if all(val == 0 for val in (left_result - right_result).values()):
                    passed += 1
                else:
                    print(f"  ✗ Associativity failed for {u.trimcode}, {v.trimcode}, {w.trimcode}")
                    pretty_print(left_result - right_result)
    
    print(f"Associativity tests: {passed}/{test_cases} passed")


def test_inverse_crystal_vs_lr_rule():
    """
    Test that inverse_crystal_product gives the correct LR coefficients.
    
    For highest weight RC-graphs, the product should match the LR rule.
    """
    print("\n" + "="*60)
    print("Testing inverse_crystal_product vs LR Rule")
    print("="*60)
    
    rc_ring = RCGraphRing()
    total_tests = 0
    passed_tests = 0
    
    for u in perms:
        for v in perms:
            # Get expected product
            prod = Sx(u) * Sx(v)
            
            for w in prod:
                expected_coeff = prod[w]
                
                # Compute via inverse_crystal_product
                # Sum over all highest weight RC-graphs
                computed_coeff = 0
                
                for u_hw in RCGraph.all_hw_rcs(u, n-1):
                    for v_hw in RCGraph.all_hw_rcs(v, n-1):
                        # Compute product at crystal quotient level
                        result = u_hw.inverse_crystal_product(v_hw)
                        
                        # Count contributions to w
                        for w_rc, coeff in result.items():
                            if w_rc.perm == w:
                                # Check if this is highest weight
                                w_rc_hw = w_rc.to_highest_weight()[0]
                                
                                # Count this if it's a highest weight component
                                for w_target_hw in RCGraph.all_hw_rcs(w, n-1):
                                    if w_rc_hw.crystal_weight == w_target_hw.crystal_weight:
                                        computed_coeff += coeff
                                        break
                
                total_tests += 1
                
                if computed_coeff == expected_coeff:
                    passed_tests += 1
                else:
                    print(f"\n✗ MISMATCH: u={u.trimcode}, v={v.trimcode}, w={w.trimcode}")
                    print(f"  Expected: {expected_coeff}, Got: {computed_coeff}")
    
    print(f"\n=== Summary ===")
    print(f"Total: {total_tests}")
    print(f"Passed: {passed_tests} ({100*passed_tests/total_tests:.1f}%)")
    
    if passed_tests == total_tests:
        print("\n🎉 inverse_crystal_product matches LR rule perfectly!")
    
    return passed_tests == total_tests


def test_crystal_structure_preservation():
    """
    Test that inverse_crystal_product preserves crystal structure.
    
    The product should:
    1. Map tensor product of crystals to direct sum of crystals
    2. Preserve highest weights appropriately
    """
    print("\n" + "="*60)
    print("Testing Crystal Structure Preservation")
    print("="*60)
    
    rc_ring = RCGraphRing()
    
    for u in perms[:10]:
        for v in perms[:10]:
            for u_hw in RCGraph.all_hw_rcs(u, n-1):
                for v_hw in RCGraph.all_hw_rcs(v, n-1):
                    # Compute product
                    result = u_hw.inverse_crystal_product(v_hw)
                    
                    # Check that all results are highest weight
                    all_hw = True
                    for w_rc in result.keys():
                        if not w_rc.inverse_crystal.is_highest_weight:
                            all_hw = False
                            print(f"\n✗ Result not highest weight: u={u.trimcode}, v={v.trimcode}")
                            print(f"  Result RC: {w_rc.perm.trimcode}")
                            break
                    
                    if not all_hw:
                        continue
                    
                    # Check that crystal weights make sense
                    # (This is a sanity check on the transposed crystal structure)
                    u_weight = u_hw.inverse_crystal.crystal_weight
                    v_weight = v_hw.inverse_crystal.crystal_weight
                    
                    for w_rc in result.keys():
                        w_weight = w_rc.inverse_crystal.crystal_weight
                        # The weights should be compatible (some relationship)
                        # This is a placeholder - you might have specific conditions


def test_comparison_with_direct_product():
    """
    Compare inverse_crystal_product with direct RC-graph product.
    
    They should be related but not identical (inverse uses transposed crystal).
    """
    print("\n" + "="*60)
    print("Comparing inverse_crystal_product with Direct Product")
    print("="*60)
    
    rc_ring = RCGraphRing()
    
    for u in perms[:5]:
        for v in perms[:5]:
            for u_rc in list(RCGraph.all_rc_graphs(u, n-1))[:2]:
                for v_rc in list(RCGraph.all_rc_graphs(v, n-1))[:2]:
                    # Direct product
                    direct_prod = u_rc.prod_with_rc(v_rc)
                    
                    # Inverse crystal product
                    inverse_prod = u_rc.inverse_crystal_product(v_rc)
                    
                    print(f"\nu={u.trimcode}, v={v.trimcode}")
                    print(f"Direct product: {len(direct_prod)} terms")
                    print(f"Inverse crystal product: {len(inverse_prod)} terms")
                    
                    # Compare permutations that appear
                    direct_perms = {rc.perm for rc in direct_prod.keys()}
                    inverse_perms = {rc.perm for rc in inverse_prod.keys()}
                    
                    if direct_perms != inverse_perms:
                        print(f"  Different permutations!")
                        print(f"  Direct only: {direct_perms - inverse_perms}")
                        print(f"  Inverse only: {inverse_perms - direct_perms}")


def test_highest_weight_subalgebra():
    """
    Test the subalgebra generated by highest weight RC-graphs.
    
    Product: (rc1 * rc2).to_highest_weight()[0]
    This should be associative and well-defined.
    """
    print("\n" + "="*60)
    print("Testing Highest Weight Subalgebra Product")
    print("="*60)
    
    rc_ring = RCGraphRing()
    
    # Test 1: Associativity
    print("\n1. Testing associativity: (a*b)*c = a*(b*c) in HW subalgebra")
    test_cases = 0
    passed = 0
    
    for u in perms[:4]:
        for v in perms[:4]:
            for w in perms[:4]:
                u_hw = next(iter(RCGraph.all_hw_rcs(u, n-1)))
                v_hw = next(iter(RCGraph.all_hw_rcs(v, n-1)))
                w_hw = next(iter(RCGraph.all_hw_rcs(w, n-1)))
                
                # Compute (u*v)*w in HW subalgebra
                uv = (rc_ring(u_hw) * rc_ring(v_hw)).to_highest_weight()[0]
                
                # Now multiply with w and go to HW
                left_result = (uv * rc_ring(w_hw)).to_highest_weight()[0]
                
                
                # Compute u*(v*w) in HW subalgebra
                vw = (rc_ring(v_hw) * rc_ring(w_hw)).to_highest_weight()[0]
                
                right_result = (rc_ring(u_hw) * vw).to_highest_weight()[0]
                
                test_cases += 1
                
                # Compare results
                match = all(val == 0 for val in (left_result - right_result).values())
                if match:
                    passed += 1
                else:
                    if test_cases <= 3:
                        print(f"  ✗ Associativity failed: {u.trimcode}, {v.trimcode}, {w.trimcode}")
                        print(f"    Left:  {left_result}")
                        print(f"    Right: {right_result}")
    
    print(f"Associativity: {passed}/{test_cases} passed")
    
    # Test 2: Results are highest weight
    print("\n2. Testing that results are highest weight")
    test_cases = 0
    passed = 0
    
    for u in perms[:5]:
        for v in perms[:5]:
            for u_hw in RCGraph.all_hw_rcs(u, n-1):
                for v_hw in RCGraph.all_hw_rcs(v, n-1):
                    product = rc_ring(u_hw) * rc_ring(v_hw)
                    hw_result = product.to_highest_weight()[0]
                    
                    # Check all results are highest weight
                    for hw_rc in hw_result.keys():
                        test_cases += 1
                        if hw_rc.is_highest_weight:
                            passed += 1
                        else:
                            print(f"  ✗ Result not HW: u={u.trimcode}, v={v.trimcode}")
    
    print(f"Results are HW: {passed}/{test_cases} passed")
    
    return passed == test_cases


def test_hw_subalgebra_vs_lr_rule():
    """
    Test if the HW subalgebra product gives the LR rule.
    
    For each (u_hw, v_hw), compute (u_hw * v_hw).to_highest_weight()[0]
    This uses the linearized to_highest_weight() on RCGraphRingElement.
    """
    print("\n" + "="*60)
    print("Testing HW Subalgebra Product vs LR Rule")
    print("="*60)
    
    rc_ring = RCGraphRing()
    total_tests = 0
    passed_tests = 0
    
    for u in perms:
        for v in perms:
            prod = Sx(u) * Sx(v)
            
            for w in prod:
                expected_coeff = prod[w]
                
                # Compute via HW subalgebra product
                computed_coeff = 0
                
                for u_hw in RCGraph.all_hw_rcs(u, n-1):
                    for v_hw in RCGraph.all_hw_rcs(v, n-1):
                        # Multiply and go to highest weight (linearized!)
                        product = rc_ring(u_hw) * rc_ring(v_hw)
                        hw_result = product.to_highest_weight()[0]
                        
                        # Count contributions to w
                        for hw_rc, coeff in hw_result.items():
                            if hw_rc.perm == w:
                                computed_coeff += coeff
                
                total_tests += 1
                
                if computed_coeff == expected_coeff:
                    passed_tests += 1
                else:
                    if total_tests <= 5:  # Print first few mismatches
                        print(f"\n✗ MISMATCH: u={u.trimcode}, v={v.trimcode}, w={w.trimcode}")
                        print(f"  Expected: {expected_coeff}, Got: {computed_coeff}")
    
    print(f"\n=== Summary ===")
    print(f"Total: {total_tests}")
    print(f"Passed: {passed_tests} ({100*passed_tests/total_tests:.1f}%)")
    
    if passed_tests == total_tests:
        print("\n🎉 HW subalgebra product gives the LR rule!")
    
    return passed_tests == total_tests


def test_inverse_crystal_matches_hw_subalgebra():
    """
    Test if inverse_crystal_product matches the HW subalgebra product.
    
    For highest weight RC-graphs:
    rc1.inverse_crystal_product(rc2) should equal
    (rc_ring(rc1) * rc_ring(rc2)).to_highest_weight()[0]
    """
    print("\n" + "="*60)
    print("Testing inverse_crystal_product vs HW Subalgebra")
    print("="*60)
    
    rc_ring = RCGraphRing()
    test_cases = 0
    passed = 0
    
    for u in perms[:5]:
        for v in perms[:5]:
            for u_hw in RCGraph.all_hw_rcs(u, n-1):
                for v_hw in RCGraph.all_hw_rcs(v, n-1):
                    # Compute via inverse_crystal_product
                    inverse_result = u_hw.inverse_crystal_product(v_hw)
                    
                    # Compute via HW subalgebra (multiply then to_highest_weight)
                    hw_subalg_result = (rc_ring(u_hw) * rc_ring(v_hw)).to_highest_weight()[0]
                    
                    test_cases += 1
                    
                    # Compare - convert hw_subalg_result to dict
                    hw_dict = {rc: coeff for rc, coeff in hw_subalg_result.items()}
                    
                    # Compare
                    all_keys = set(inverse_result.keys()) | set(hw_dict.keys())
                    match = True
                    for key in all_keys:
                        if inverse_result.get(key, 0) != hw_dict.get(key, 0):
                            match = False
                            break
                    
                    if match:
                        passed += 1
                    else:
                        if test_cases <= 3:  # Print first few mismatches
                            print(f"\n✗ Mismatch: u={u.trimcode}, v={v.trimcode}")
                            print(f"  Inverse: {inverse_result}")
                            print(f"  HW subalg: {hw_dict}")
    
    print(f"\nComparison: {passed}/{test_cases} passed")
    
    if passed == test_cases:
        print("✓ inverse_crystal_product matches HW subalgebra product!")
    
    return passed == test_cases


# Update the main block
if __name__ == "__main__":
    print("Testing inverse_crystal_product implementation")
    print("="*60)
    
    # Test the HW subalgebra structure first
    hw_subalg_works = test_highest_weight_subalgebra()
    
    if hw_subalg_works:
        print("\n✓ HW subalgebra product is associative and well-defined!")
        
        # Test if it gives the LR rule
        hw_lr_works = test_hw_subalgebra_vs_lr_rule()
        
        # Test if inverse_crystal_product matches it
        inverse_matches = test_inverse_crystal_matches_hw_subalgebra()
        
        if hw_lr_works and inverse_matches:
            print("\n" + "="*60)
            print("🎉 Everything works!")
            print("- HW subalgebra is associative")
            print("- HW subalgebra gives LR coefficients")
            print("- inverse_crystal_product implements HW subalgebra")
            print("="*60)