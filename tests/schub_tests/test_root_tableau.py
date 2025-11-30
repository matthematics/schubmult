"""
Test suite comparing direct crystal operators against the round-trip RC-graph operators.
"""
import pytest
from schubmult import Permutation, RCGraph, RootTableau


def test_raising_operator_direct_vs_roundtrip():
    """Compare raising_operator_direct against raising_operator for several cases."""
    test_cases = []
    
    # Generate test cases
    for n in range(3, 6):
        for perm in Permutation.all_permutations(n):
            rc = RCGraph.principal_rc(perm, n - 1)
            rt = RootTableau.from_rc_graph(rc)
            test_cases.append((rt, perm))
    
    failures = []
    successes = 0
    
    for rt, perm in test_cases:
        for i in range(1, len(perm.trimcode)):
            # Try both operators (assuming you have raising_operator implemented)
            result_roundtrip = rt.raising_operator(i)
            
            result_direct = rt.raising_operator_direct(i)
            
            # Both should be None or both should succeed
            if result_roundtrip is None and result_direct is None:
                successes += 1
                continue
            
            if result_roundtrip is None or result_direct is None:
                failures.append({
                    'perm': perm,
                    'i': i,
                    'roundtrip': result_roundtrip,
                    'direct': result_direct,
                    'reason': f'result_roundtrip returned {result_roundtrip}, other returned {result_direct}'
                })
                continue
            
            # Compare the results
            if result_roundtrip != result_direct:
                failures.append({
                    'perm': perm,
                    'i': i,
                    'roundtrip': result_roundtrip,
                    'direct': result_direct,
                    'reason': 'Results differ'
                })
            else:
                successes += 1
    
    if failures:
        print(f"\n{len(failures)} failures out of {len(failures) + successes} tests:")
        for failure in failures[:5]:
            print(f"  perm={failure['perm']}, i={failure['i']}: {failure['reason']}")
    
    assert len(failures) == 0, f"{len(failures)} failures (see output)"


def test_lowering_operator_direct_vs_roundtrip():
    """Compare lowering_operator_direct against lowering_operator for several cases."""
    test_cases = []
    
    # Generate test cases
    for n in range(3, 6):
        for perm in Permutation.all_permutations(n):
            rc = RCGraph.principal_rc(perm, n - 1).to_highest_weight()[0] # needs to be highest weight because principal is already lowest weight
            rt = RootTableau.from_rc_graph(rc)
            test_cases.append((rt, perm))
    
    failures = []
    successes = 0
    
    for rt, perm in test_cases:
        for i in range(1, len(perm.trimcode)):
            # Try both operators (assuming you have lowering_operator implemented)
            result_roundtrip = rt.lowering_operator(i)
            
            result_direct = rt.lowering_operator_direct(i)
            
            # Both should be None or both should succeed
            if result_roundtrip is None and result_direct is None:
                successes += 1
                continue
            
            if result_roundtrip is None or result_direct is None:
                failures.append({
                    'perm': perm,
                    'i': i,
                    'roundtrip': result_roundtrip,
                    'direct': result_direct,
                    'reason': f'result_roundtrip returned {result_roundtrip}, other returned {result_direct}'
                })
                continue
            
            # Compare the results
            if result_roundtrip != result_direct:
                failures.append({
                    'perm': perm,
                    'i': i,
                    'roundtrip': result_roundtrip,
                    'direct': result_direct,
                    'reason': 'Results differ'
                })
            else:
                successes += 1
    
    if failures:
        print(f"\n{len(failures)} failures out of {len(failures) + successes} tests:")
        for failure in failures[:5]:
            print(f"  perm={failure['perm']}, i={failure['i']}: {failure['reason']}")
    
    assert len(failures) == 0, f"{len(failures)} failures (see output)"


def test_direct_operator_preserves_invariants():
    """Verify that direct operators preserve the key invariants."""
    for n in range(3, 5):
        for perm in Permutation.all_permutations(n):
            rc = RCGraph.principal_rc(perm, n - 1)
            rt = RootTableau.from_rc_graph(rc)
            
            for i in range(1, len(perm.trimcode)):
                result = rt.raising_operator_direct(i)
                if result is None:
                    continue
                
                # Check invariants
                assert result.perm == rt.perm, f"Permutation changed: {rt.perm} -> {result.perm}"
                assert result.shape == rt.shape, f"Shape changed: {rt.shape} -> {result.shape}"
                #assert result.rc_graph == rt.rc_graph, f"RC-graph changed"
                
                # Check EG-invariant
                assert result.edelman_greene_invariant == rt.edelman_greene_invariant, \
                    f"EG-invariant changed: {rt.edelman_greene_invariant} -> {result.edelman_greene_invariant}"

def test_original_operator_preserves_invariants():
    """Verify that direct operators preserve the key invariants."""
    for n in range(3, 5):
        for perm in Permutation.all_permutations(n):
            rc = RCGraph.principal_rc(perm, n - 1)
            rt = RootTableau.from_rc_graph(rc)
            
            for i in range(1, len(perm.trimcode)):
                result = rt.raising_operator(i)
                if result is None:
                    continue
                
                # Check invariants
                assert result.perm == rt.perm, f"Permutation changed: {rt.perm} -> {result.perm}"
                assert result.shape == rt.shape, f"Shape changed: {rt.shape} -> {result.shape}"
                #assert result.rc_graph == rt.rc_graph, f"RC-graph changed"
                
                # Check EG-invariant
                assert result.edelman_greene_invariant == rt.edelman_greene_invariant, \
                    f"EG-invariant changed: {rt.edelman_greene_invariant} -> {result.edelman_greene_invariant}"



def test_direct_operator_simple_case():
    """Test a single simple case manually for debugging."""
    perm = Permutation([2, 1, 3])
    rc = RCGraph.principal_rc(perm, 2)
    rt = RootTableau.from_rc_graph(rc)
    
    print("\nOriginal root tableau:")
    print(rt)
    
    for i in range(1, len(perm.trimcode)):
        result_direct = rt.raising_operator_direct(i)
        result_roundtrip = rt.raising_operator(i)
        
        print(f"\ne_{i} direct:")
        print(result_direct)
        print(f"\ne_{i} roundtrip:")
        print(result_roundtrip)
        
        if result_direct is not None and result_roundtrip is not None:
            assert result_direct == result_roundtrip, f"Mismatch for e_{i}"


if __name__ == "__main__":
    # Run a simple test for quick debugging
    test_direct_operator_simple_case()
    print("\n" + "="*60)
    print("Running full test suite...")
    pytest.main([__file__, "-v"])