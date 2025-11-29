from sympy import init_printing, pretty_print
import sympy

from schubmult.schub_lib.perm_lib import Permutation
from schubmult.schub_lib.rc_graph import RCGraph
from schubmult.schub_lib.crystal_graph import CrystalGraphTensor
from schubmult.rings.schubert_ring import Sx
from schubmult.rings.rc_graph_ring import RCGraphRing

init_printing()

# n = 3  # Rank

# Get all permutations


def dominifiable(u, v):
    v_work = v
    while v_work != v.minimal_dominant_above():
        index = max([i for i in range(len(v_work.trimcode)-1) if v_work.trimcode[i] < v_work.trimcode[i+1]])
        if index in u.descents():
            return False
        v_work = v_work.swap(index, index+1)
    return True

def is_isomorphic_crystal(tensor, rc_graph):
    hw1 = tensor.to_highest_weight()[0]
    hw2 = rc_graph.to_highest_weight()[0]
    if hw1.crystal_weight != hw2.crystal_weight:
        return False
    lw1 = hw1.to_lowest_weight()[0]
    lw2 = hw2.to_lowest_weight()[0]
    if lw1.crystal_weight != lw2.crystal_weight:
        return False
    return True


def tensor_decomposes(u_rc0, v_rc0):
    """
    Returns True if Crystal(u_rc0) ⊗ Crystal(v_rc0) decomposes into 
    a direct sum of Demazure crystals.
    
    Method: Check if every tensor element is reachable from some highest weight.
    """
    # Collect all tensor elements
    all_tensors = set()
    for u_rc in u_rc0.full_crystal:
        for v_rc in v_rc0.full_crystal:
            all_tensors.add(CrystalGraphTensor(u_rc, v_rc))
    
    # Find all highest weight components
    visited = set()
    components = []
    
    for tensor in all_tensors:
        if tensor in visited:
            continue
        
        # Go to highest weight
        hw = tensor.to_highest_weight()[0]
        
        # Get full crystal component from this highest weight
        component = hw.full_crystal
        components.append(hw)
        visited.update(component)
    
    # Decomposes ⟺ we visited everything
    return len(visited) == len(all_tensors)


def is_demazure_crystal(component):
    """
    Check if a crystal component is actually a Demazure crystal.
    
    A Demazure crystal must have:
    - Unique highest weight
    - Unique lowest weight
    """
    highest_weights = set()
    lowest_weights = set()
    
    for element in component:
        # Check if this is a highest weight (all e_i return None)
        is_hw = True
        for i in range(1, element.crystal_length()):
            if element.raising_operator(i) is not None:
                is_hw = False
                break
        if is_hw:
            highest_weights.add(element)
        
        # Check if this is a lowest weight (all f_i return None)
        is_lw = True
        for i in range(1, element.crystal_length()):
            if element.lowering_operator(i) is not None:
                is_lw = False
                break
        if is_lw:
            lowest_weights.add(element)
    
    # Demazure crystal has exactly one highest and one lowest weight
    return len(highest_weights) == 1 and len(lowest_weights) == 1


def tensor_decomposes_into_demazure(u_rc0, v_rc0):
    """
    Returns True if Crystal(u_rc0) ⊗ Crystal(v_rc0) decomposes into 
    a direct sum of DEMAZURE crystals.
    
    Checks:
    1. All elements are reachable (components partition the tensor)
    2. Each component is a Demazure crystal (unique highest/lowest weight)
    """
    # Collect all tensor elements
    all_tensors = set()
    for u_rc in u_rc0.full_crystal:
        for v_rc in v_rc0.full_crystal:
            all_tensors.add(CrystalGraphTensor(u_rc, v_rc))
    
    # Find all components
    visited = set()
    components = []
    
    for tensor in all_tensors:
        if tensor in visited:
            continue
        
        # Go to highest weight
        hw = tensor.to_highest_weight()[0]
        
        # Get full crystal component from this highest weight
        component = set(hw.full_crystal)
        components.append(component)
        visited.update(component)
    
    # Check 1: All elements covered (already known to pass)
    if len(visited) != len(all_tensors):
        return False
    
    # Check 2: Each component is a Demazure crystal
    for component in components:
        if not is_demazure_crystal(component):
            # print(f"  Component is NOT Demazure (size {len(component)})")
            # Debug: how many highest/lowest weights?
            hw_count = sum(1 for e in component if all(
                e.raising_operator(i) is None 
                for i in range(1, e.crystal_length())
            ))
            lw_count = sum(1 for e in component if all(
                e.lowering_operator(i) is None 
                for i in range(1, e.crystal_length())
            ))
            # print(f"    Highest weights: {hw_count}")
            # print(f"    Lowest weights: {lw_count}")
            return False
    
    return True


def test_decomposition(dominifiable_case=True):
    """
    Test whether tensor products decompose into Demazure crystals
    """
    decompose_count = 0
    dont_decompose_count = 0
    non_demazure_examples = []
    
    for u in perms:
        for v in perms:
            if dominifiable_case and not dominifiable(u, v):
                continue
            
            for u_rc in RCGraph.all_hw_rcs(u, n-1):
                for v_rc in RCGraph.all_hw_rcs(v, n-1):
                    decomposes = tensor_decomposes_into_demazure(u_rc, v_rc)
                    
                    if decomposes:
                        decompose_count += 1
                    else:
                        dont_decompose_count += 1
                        non_demazure_examples.append({
                            'u': u.trimcode,
                            'v': v.trimcode,
                            'u_rc': u_rc,
                            'v_rc': v_rc
                        })
                        # print(f"\n!!! NON-DEMAZURE COMPONENT FOUND !!!")
                        # print(f"u={u.trimcode}, v={v.trimcode}")
    
    print(f"\n=== Summary ===")
    print(f"Decomposes into Demazure: {decompose_count}")
    print(f"Has non-Demazure components: {dont_decompose_count}")
    print(f"Dominifiable={dominifiable_case}")
    if len(non_demazure_examples) > 0:
        print(f"\nFound {len(non_demazure_examples)} cases with non-Demazure components!")
    else:
        print("\nAll tensor products decompose into Demazure crystals!")
    
    return non_demazure_examples


def interleave_rc_graphs(u_rc, v_rc, pattern='alternating'):
    """
    Interleave two RC-graphs row by row.
    
    Args:
        u_rc, v_rc: RC graphs to interleave
        pattern: 'alternating' or 'u_first' (concatenation)
    
    Returns:
        List of rows to be reduced via RC-graph product
    """
    if pattern == 'alternating':
        # Alternate rows: u1, v1, u2, v2, ...
        result_rows = []
        max_len = max(len(u_rc), len(v_rc))
        
        for i in range(max_len):
            if i < len(u_rc):
                result_rows.append(u_rc[i])
            if i < len(v_rc):
                result_rows.append(v_rc[i])
        
        return result_rows
    
    elif pattern == 'u_first':
        # Simple concatenation: all u, then all v
        return list(u_rc) + list(v_rc)
    
    else:
        raise NotImplementedError(f"Pattern {pattern} not implemented")


def reduce_interleaved(interleaved_rows, n):
    """
    Reduce an interleaved list of rows via iterated RC-graph product.
    
    Start with first row, multiply with second, then third, etc.
    Returns set of all resulting RC-graphs.
    """
    if len(interleaved_rows) == 0:
        return {RCGraph([])}
    
    # Start with the first row as an RC-graph
    current_rcs = {RCGraph([interleaved_rows[0]])}
    
    # Multiply row by row
    for row in interleaved_rows[1:]:
        next_rcs = set()
        row_rc = RCGraph([row])
        
        for curr_rc in current_rcs:
            # Multiply curr_rc with row_rc
            product = curr_rc.prod_with_rc(row_rc)
            next_rcs.update(product)
        
        current_rcs = next_rcs
    
    return current_rcs


def test_interleaving_lr_rule(pattern='alternating'):
    """
    Test if interleaving + reduction gives the LR rule.
    
    Pattern can be 'alternating' or 'u_first' (simple concatenation).
    """
    print("\n" + "="*60)
    print(f"Testing Interleaving LR Rule (pattern={pattern})")
    print("="*60)
    
    total_tests = 0
    passed_tests = 0
    failed_tests = 0
    
    for u in perms:
        for v in perms:
            prod = Sx(u) * Sx(v)
            
            for w in prod:
                expected_coeff = prod[w]
                
                # Count via interleaving rule
                computed_coeff = 0
                
                for u_rc in RCGraph.all_rc_graphs(u, n-1):
                    for v_rc in RCGraph.all_rc_graphs(v, n-1):
                        # Interleave
                        interleaved = interleave_rc_graphs(u_rc, v_rc, pattern=pattern)
                        
                        # Reduce via RC-graph product
                        result_rcs = reduce_interleaved(interleaved, n)
                        
                        # Check if w appears in the result
                        for result_rc in result_rcs:
                            if result_rc.perm == w:
                                computed_coeff += 1
                
                total_tests += 1
                
                # Check if correct
                if computed_coeff == expected_coeff:
                    passed_tests += 1
                else:
                    failed_tests += 1
                    print(f"\n✗ MISMATCH: u={u.trimcode}, v={v.trimcode}, w={w.trimcode}")
                    print(f"  Expected: {expected_coeff}, Got: {computed_coeff}")
                    
                    if computed_coeff > expected_coeff:
                        print(f"  Over-counted by {computed_coeff - expected_coeff}")
                    else:
                        print(f"  Under-counted by {expected_coeff - computed_coeff}")
    
    print(f"\n=== Summary for pattern={pattern} ===")
    print(f"Total tests: {total_tests}")
    print(f"Passed: {passed_tests} ({100*passed_tests/total_tests:.1f}%)")
    print(f"Failed: {failed_tests} ({100*failed_tests/total_tests:.1f}%)")
    
    if failed_tests == 0:
        print("\n🎉 SUPER ASSOCIATIVITY CONFIRMED! 🎉")
        print("Interleaving gives a fully general bijective LR rule!")
    else:
        print(f"\nSuper associativity does not hold for pattern={pattern}")
        print("But we've learned something about where it fails!")
    
    return passed_tests == total_tests


def transposed_crystal_operators(rc, i):
    """
    Apply crystal operators using the TRANSPOSED crystal structure.
    
    Transpose -> apply operator -> transpose back to original dimensions.
    
    Args:
        rc: RC-graph
        i: crystal operator index
    
    Returns:
        (raising_result, lowering_result) where each is the RC-graph after
        applying the transposed crystal operator (or None if not applicable)
    """
    original_len = len(rc)
    
    # Transpose (length doesn't matter for first transpose)
    rc_t = rc.transpose()
    
    # Apply raising operator in transposed space
    raised_t = rc_t.raising_operator(i)
    raised = raised_t.transpose(original_len) if raised_t is not None else None
    
    # Apply lowering operator in transposed space
    lowered_t = rc_t.lowering_operator(i)
    lowered = lowered_t.transpose(original_len) if lowered_t is not None else None
    
    return raised, lowered


def tensor_decomposes_transposed(u_rc0, v_rc0):
    """
    Check if Crystal(u_rc0) ⊗ Crystal(v_rc0) decomposes when using
    TRANSPOSED crystal structure on the product range.
    
    The RC-graph product is a crystal homomorphism when:
    1. Domain has tensor product crystal structure
    2. Range has TRANSPOSED crystal structure
    """
    # Compute all products (with transposed crystal structure)
    ring = RCGraphRing()
    all_products_t = set()
    
    for u_rc in u_rc0.full_crystal:
        for v_rc in v_rc0.full_crystal:
            product = ring(u_rc) * ring(v_rc)
            
            # Add transposed versions to our set
            for result_rc in product.keys():
                # Transpose to get the crystal structure we want
                result_t = result_rc.transpose()
                all_products_t.add(result_t)
    
    # Find connected components in TRANSPOSED crystal space
    visited = set()
    components = []
    
    for prod_t in all_products_t:
        if prod_t in visited:
            continue
        
        # Find highest weight in TRANSPOSED space
        hw_t = prod_t
        changed = True
        while changed:
            changed = False
            for i in range(1, prod_t.crystal_length()):
                raised = hw_t.raising_operator(i)
                if raised is not None:
                    hw_t = raised
                    changed = True
                    break
        
        # Get full crystal component from this highest weight (in transposed space)
        component_t = set()
        to_explore = [hw_t]
        explored = set()
        
        while to_explore:
            current = to_explore.pop()
            if id(current) in explored:
                continue
            explored.add(id(current))
            component_t.add(current)
            
            # Apply all lowering operators in transposed space
            for i in range(1, current.crystal_length()):
                lowered = current.lowering_operator(i)
                if lowered is not None and id(lowered) not in explored:
                    to_explore.append(lowered)
        
        components.append(component_t)
        visited.update(component_t)
    
    # Check if we covered everything
    return len(visited) == len(all_products_t)


def is_demazure_crystal_transposed(component_t):
    """
    Check if a component in TRANSPOSED crystal space is a Demazure crystal.
    
    Args:
        component_t: Set of RC-graphs in transposed crystal space
    
    Returns:
        True if unique highest and lowest weight in transposed space
    """
    highest_weights = set()
    lowest_weights = set()
    
    for element in component_t:
        # Check if highest weight in transposed space
        is_hw = True
        for i in range(1, element.crystal_length()):
            if element.raising_operator(i) is not None:
                is_hw = False
                break
        if is_hw:
            highest_weights.add(element)
        
        # Check if lowest weight in transposed space
        is_lw = True
        for i in range(1, element.crystal_length()):
            if element.lowering_operator(i) is not None:
                is_lw = False
                break
        if is_lw:
            lowest_weights.add(element)
    
    return len(highest_weights) == 1 and len(lowest_weights) == 1


def tensor_decomposes_into_demazure_transposed(u_rc0, v_rc0):
    """
    Check if tensor decomposes into Demazure crystals using TRANSPOSED structure.
    """
    # Compute all products in transposed space
    ring = RCGraphRing()
    all_products_t = set()
    
    for u_rc in u_rc0.full_crystal:
        for v_rc in v_rc0.full_crystal:
            product = ring(u_rc) * ring(v_rc)
            for result_rc in product.keys():
                all_products_t.add(result_rc.transpose())
    
    # Find components in transposed space
    visited = set()
    components = []
    
    for prod_t in all_products_t:
        if prod_t in visited:
            continue
        
        # Get highest weight in transposed space
        hw_t = prod_t
        changed = True
        while changed:
            changed = False
            for i in range(1, prod_t.crystal_length()):
                raised = hw_t.raising_operator(i)
                if raised is not None:
                    hw_t = raised
                    changed = True
                    break
        
        # Build component from highest weight
        component_t = set()
        to_explore = [hw_t]
        explored = set()
        
        while to_explore:
            current = to_explore.pop()
            if id(current) in explored:
                continue
            explored.add(id(current))
            component_t.add(current)
            
            for i in range(1, current.crystal_length()):
                lowered = current.lowering_operator(i)
                if lowered is not None and id(lowered) not in explored:
                    to_explore.append(lowered)
        
        components.append(component_t)
        visited.update(component_t)
    
    # Check coverage
    if len(visited) != len(all_products_t):
        return False
    
    # Check if each component is Demazure
    for component_t in components:
        if not is_demazure_crystal_transposed(component_t):
            return False
    
    return True


def test_transposed_decomposition(dominifiable_case=True):
    """
    Test whether tensor products decompose into Demazure crystals
    when using TRANSPOSED crystal structure.
    """
    print("\n" + "="*60)
    print("Testing Decomposition with TRANSPOSED Crystal Structure")
    print(f"Dominifiable cases only: {dominifiable_case}")
    print("="*60)
    
    decompose_count = 0
    dont_decompose_count = 0
    
    for u in perms:
        for v in perms:
            if dominifiable_case and not dominifiable(u, v):
                continue
            
            for u_rc in RCGraph.all_hw_rcs(u, n-1):
                for v_rc in RCGraph.all_hw_rcs(v, n-1):
                    decomposes = tensor_decomposes_into_demazure_transposed(u_rc, v_rc)
                    
                    if decomposes:
                        decompose_count += 1
                    else:
                        dont_decompose_count += 1
                        print(f"\n✗ Non-Demazure (transposed): u={u.trimcode}, v={v.trimcode}")
    
    print(f"\n=== Summary (Transposed) ===")
    print(f"Decomposes into Demazure: {decompose_count}")
    print(f"Non-Demazure components: {dont_decompose_count}")
    
    if dont_decompose_count == 0:
        print("\n🎉 TRANSPOSED crystal structure gives Demazure decomposition!")
        print("This means the RC-graph product is a crystal homomorphism!")
    
    return dont_decompose_count == 0


def test_divided_difference_preserves_crystal():
    """
    Test that divided difference preserves crystal structure in dominifiable case.
    
    The divided difference should:
    1. Not affect the u part (by construction)
    2. Map crystals bijectively: w_dom_crystal → w_crystal
    3. Preserve crystal weights
    """
    print("\n" + "="*60)
    print("Testing Divided Difference Preserves Crystal Structure")
    print("="*60)
    
    for u in perms:
        for v in perms:
            if not dominifiable(u, v):
                continue
            
            dom_perm = v.minimal_dominant_above()
            diff_elem = (~v) * dom_perm
            
            # For each w_dom in u * dom_perm
            prod_dom = Sx(u) * Sx(dom_perm)
            
            for w_dom in prod_dom:
                w = w_dom * (~diff_elem)
                
                if w.inv != w_dom.inv - diff_elem.inv:
                    continue
                
                # For each w_dom_rc crystal
                for w_dom_rc in RCGraph.all_hw_rcs(w_dom, n-1):
                    w_dom_weight = w_dom_rc.crystal_weight
                    w_dom_crystal_size = len(list(w_dom_rc.full_crystal))
                    
                    # Apply divided difference to get w_rcs
                    w_rcs = w_dom_rc.divdiff_perm(diff_elem)
                    
                    if len(w_rcs) == 0:
                        print(f"WARNING: divdiff kills w_dom_rc for w_dom={w_dom.trimcode}")
                        continue
                    
                    print(f"\nw_dom={w_dom.trimcode} -> w={w.trimcode}")
                    print(f"  divdiff produces {len(w_rcs)} w_rcs")
                    
                    # Check crystal sizes match
                    for w_rc in w_rcs:
                        w_crystal_size = len(list(w_rc.full_crystal))
                        w_rc_weight = w_rc.crystal_weight
                        
                        print(f"    w_rc: weight={w_rc_weight}, crystal_size={w_crystal_size}")
                        
                        if w_crystal_size == w_dom_crystal_size:
                            print(f"      ✓ Crystal sizes match: {w_crystal_size}")
                            
                            # Check if weights match
                            if w_rc_weight == w_dom_weight:
                                print(f"      ✓ Crystal weights match: {w_dom_weight}")
                                
                                # Check if they're isomorphic
                                if is_isomorphic_crystal(w_dom_rc, w_rc):
                                    print(f"      ✓ Crystals are isomorphic!")
                        else:
                            print(f"      ✗ Crystal size mismatch: {w_dom_crystal_size} -> {w_crystal_size}")


def test_crystal_correspondence_dominifiable():
    """
    Test the crystal-level correspondence in the dominifiable case.
    
    For u * v where v is dominifiable:
    1. Raise v to dom_perm
    2. Crystal(u_rc) ⊗ Crystal(dom_rc) decomposes
    3. Each component is isomorphic to a w_rc crystal
    4. Apply divided difference: should preserve the correspondence
    """
    print("\n" + "="*60)
    print("Testing Crystal Correspondence in Dominifiable Case")
    print("="*60)
    
    for u in perms:
        for v in perms:
            if not dominifiable(u, v):
                continue
            
            dom_perm = v.minimal_dominant_above()
            diff_elem = (~v) * dom_perm
            
            print(f"\nu={u.trimcode}, v={v.trimcode}, dom={dom_perm.trimcode}")
            print(f"diff_elem={diff_elem.trimcode} (length {diff_elem.inv})")
            
            # Compute both products
            prod_dom = Sx(u) * Sx(dom_perm)  # Dominant case
            prod_v = Sx(u) * Sx(v)           # Target product
            
            for w_dom in prod_dom:
                # Apply divided difference to get w
                w = w_dom * (~diff_elem)
                
                if w.inv != w_dom.inv - diff_elem.inv:
                    # Divided difference kills this term
                    continue
                
                print(f"\n  w_dom={w_dom.trimcode} -> w={w.trimcode}")
                
                # For each u_rc that contributes to w_dom
                for u_rc in RCGraph.all_hw_rcs(u, n-1):
                    for dom_rc in RCGraph.all_hw_rcs(dom_perm, n-1):
                        # Check if this (u_rc, dom_rc) contributes to w_dom
                        dp_ret = u_rc.dualpieri(dom_perm, w_dom)
                        
                        if len(dp_ret) == 0:
                            continue
                        
                        print(f"    Found contribution from u_rc, dom_rc via dualpieri")
                        
                        # For each w_dom_rc that dualpieri found
                        for coeff, w_dom_rc in dp_ret:
                            print(f"      w_dom_rc: weight={w_dom_rc.crystal_weight}")
                            
                            # Apply divided difference to w_dom_rc
                            w_rcs = w_dom_rc.divdiff_perm(diff_elem)
                            
                            print(f"        divdiff produces {len(w_rcs)} w_rcs")
                            
                            for w_rc in w_rcs:
                                print(f"          w_rc: perm={w_rc.perm.trimcode}, weight={w_rc.crystal_weight}")
                                
                                # This w_rc should contribute to c^w_{u,v}
                                if w_rc.perm == w:
                                    print(f"          ✓ Correct permutation!")
                                    
                                    # Check if it's a valid contribution
                                    # (Should we verify it appears in the actual product?)
                                    expected_in_prod = w in prod_v
                                    if expected_in_prod:
                                        print(f"          ✓ w is in u*v product")
                                else:
                                    print(f"          ✗ Wrong permutation: got {w_rc.perm.trimcode}, expected {w.trimcode}")


def verify_dominifiable_lr_rule_crystal_level():
    """
    Verify the complete dominifiable LR rule at the crystal level.
    
    The rule should be:
    c^w_{u,v} = number of contributions where:
    1. (u_rc, dom_rc) contributes to w_dom via dualpieri
    2. w_dom_rc.divdiff_perm(diff_elem) gives w_rcs
    3. All w_rcs from same w_dom_rc should have same highest weight
    4. Count one contribution per (u_rc, dom_rc, w_hw) triple
    """
    print("\n" + "="*60)
    print("Verifying Dominifiable LR Rule at Crystal Level")
    print("="*60)
    
    for u in perms:
        for v in perms:
            if not dominifiable(u, v):
                continue
            
            dom_perm = v.minimal_dominant_above()
            diff_elem = (~v) * dom_perm
            prod = Sx(u) * Sx(v)
            
            for w in prod:
                expected_coeff = prod[w]
                
                # Find w_dom
                w_dom = w * diff_elem
                if Sx(u) * Sx(dom_perm) == 0 or w_dom not in (Sx(u) * Sx(dom_perm)):
                    continue
                
                # Count at crystal level using divdiff_perm
                computed_coeff = 0
                
                for u_rc in RCGraph.all_hw_rcs(u, n-1):
                    for dom_rc in RCGraph.all_hw_rcs(dom_perm, n-1):
                        # Check if contributes to w_dom
                        dp_ret = u_rc.dualpieri(dom_perm, w_dom)
                        
                        for coeff_dp, w_dom_rc in dp_ret:
                            # Apply divided difference
                            w_rcs = w_dom_rc.divdiff_perm(diff_elem)
                            
                            # Check that all w_rcs with perm == w have same highest weight
                            w_matching_rcs = [rc for rc in w_rcs if rc.perm == w]
                            
                            if len(w_matching_rcs) > 0:
                                # Get highest weights
                                hw_weights = set()
                                for w_rc in w_matching_rcs:
                                    w_hw = w_rc.to_highest_weight()[0]
                                    hw_weights.add(w_hw.crystal_weight)
                                
                                if len(hw_weights) > 1:
                                    print(f"WARNING: Multiple highest weights for w={w.trimcode}")
                                    print(f"  u={u.trimcode}, v={v.trimcode}")
                                    print(f"  Highest weights: {hw_weights}")
                                
                                # Count this as ONE contribution (one crystal component)
                                computed_coeff += 1
                
                if computed_coeff == expected_coeff:
                    print(f"✓ u={u.trimcode}, v={v.trimcode}, w={w.trimcode}: {computed_coeff}")
                else:
                    print(f"✗ MISMATCH: u={u.trimcode}, v={v.trimcode}, w={w.trimcode}")
                    print(f"  Expected: {expected_coeff}, Got: {computed_coeff}")
                    if computed_coeff > expected_coeff:
                        print(f"  Over-counted by {computed_coeff - expected_coeff}")
                    else:
                        print(f"  Under-counted by {expected_coeff - computed_coeff}")


# Add to main block
if __name__ == "__main__":
    # print("Testing Dominifiable LR Rule")
    # print("="*60)
    
    # # Test 1: Check decomposition
    # # print("\n1. Testing tensor decomposition (dominifiable cases):")
    # # non_demazure_examples = test_decomposition(dominifiable_case=True)
    
    # # print("\n2. Testing tensor decomposition (all cases):")
    # # all_non_demazure = test_decomposition(dominifiable_case=False)
    
    # # # Test 2: Interleaving
    # # print("\n3. Testing interleaving (u_first):")
    # # test_interleaving_lr_rule(pattern='u_first')
    
    # # print("\n4. Testing interleaving (alternating):")
    # # test_interleaving_lr_rule(pattern='alternating')
    
    # # # Test 3: Transposed crystal
    # # print("\n5. Testing transposed crystal (dominifiable):")
    # # test_transposed_decomposition(dominifiable_case=True)
    
    # # print("\n6. Testing transposed crystal (all):")
    # # test_transposed_decomposition(dominifiable_case=False)
    
    # # Test 4: Crystal correspondence
    # print("\n7. Testing crystal correspondence:")
    # test_crystal_correspondence_dominifiable()
    
    # print("\n8. Testing divided difference preservation:")
    # test_divided_difference_preserves_crystal()
    
    # print("\n9. Verifying LR rule at crystal level:")
    # verify_dominifiable_lr_rule_crystal_level()

    # middle insertion
    # import sys
    # n = int(sys.argv[1])
    # perms = Permutation.all_permutations(n)
    # for w in perms:
    #     # for rc in RCGraph.all_rc_graphs(w, n - 1):
    #     #     print("w_rc")
    #     #     pretty_print(rc)
    #     #     start = 1
    #     #     end  = len(rc) - 1
    #     #     rc_ring = RCGraphRing()
    #     #     middle_rc, bottom_rc = rc.rowrange(start).vertical_cut(end - start)
    #     #     top_rc = rc.vertical_cut(start)[0]
    #     #     print(Sx(w).coproduct(*list(range(1, n - 1))))
    #     #     print("top * bottom")
    #     #     r_elem1 = rc_ring(top_rc) * rc_ring(bottom_rc)
    #     #     pretty_print(r_elem1)
    #     #     pretty_print(middle_rc)
    #     #     r_elem = rc_ring(top_rc) * rc_ring(middle_rc) * rc_ring(bottom_rc)
    #     #     print("top * middle * bottom coeff")
    #     #     print(r_elem.get(rc))
    #     #     print("Coeff in coprod")
    #     tm = len(w.trimcode)
    #     if tm < 4:
    #         continue
    #     copfot = Sx(w).coproduct(*list(range(2, tm)))
    #     fopcot = {}
    #     print(f"{copfot}")
        rc_ring = RCGraphRing()

    #     for (perm11, perm22), coeff in copfot.items():
    #         perm1, perm2 = Permutation(perm11), Permutation(perm22)
    #         working_perm2 = perm2
    #         working_w = w
    #         while len(working_perm2.trimcode) > 1:
    #             d = max(working_perm2.descents())
    #             working_perm2 = working_perm2.swap(d, d + 1)
    #             working_w = working_w.swap(d, d + 1)
    #         rc1 = RCGraph.principal_rc(perm1, tm - 2)
    #         rc2 = RCGraph.principal_rc(working_perm2, 1)
    #         for w_rc in RCGraph.all_rc_graphs(working_w, tm):
    #             if w_rc.vertical_cut(1) == (rc2, rc1.extend(1)):
    #                 fopcot[(perm1, perm2)] = fopcot.get((perm1, perm2), 0) + 1
    #     print(f"fopcot: {fopcot}")
    #     for key in copfot:
    #         if copfot[key] != fopcot.get(key, 0):
    #             print(f"Mismatch for {key}: copfot={copfot[key]}, fopcot={fopcot.get((Permutation(key[0]), Permutation(key[1])), 0)}")
    #         else:
    #             print(f"Match for {key}: {copfot[key]}")
            
            