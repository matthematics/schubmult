from sympy import init_printing, pretty_print
import sympy

def dominifiable(u, v):
    v_work = v
    while v_work != v.minimal_dominant_above():
        index = max([i for i in range(len(v_work.trimcode)-1) if v_work.trimcode[i] < v_work.trimcode[i+1]])
        if index in u.descents():
            return False
        v_work = v_work.swap(index, index+1)
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

def crystal_dom_product(dom_rc, u):
    from schubmult import RCGraphRing, Sx
    ring = RCGraphRing()
    res = ring.zero
    cheat_prod = Sx(dom_rc.perm) * Sx(u)
    # if len(cheat_prod) == 1:
    #     w = next(iter(cheat_prod))
    #     pants =  CrystalGraphTensor(dom_rc, u_rc_lw).to_lowest_weight()[0]
    #     wt = pants.to_lowest_weight()[0].crystal_weight
    #     w_prin = RCGraph.principal_rc(w, len(dom_rc))
    #     if w_prin.crystal_weight == wt:
    #         return ring(w_prin)
    #     return res
    for w in cheat_prod:
        for u_rc_hw in RCGraph.all_hw_rcs(u, n-1):
            for u_rc in u_rc_hw.full_crystal:
                    dp_ret = u_rc.dualpieri(dom_rc.perm, w)
                    if len(dp_ret) > 0:
                        pants =  CrystalGraphTensor(dom_rc, u_rc).to_lowest_weight()[0]
                        wt = pants.to_lowest_weight()[0].crystal_weight
                        wp_rcs = [rc for rc in RCGraph.all_rc_graphs(w, n-1, weight=wt) if rc.is_lowest_weight]
                        wp_rc = wp_rcs[0]
                    
                        if wp_rc.to_highest_weight()[0].crystal_weight == pants.to_highest_weight()[0].crystal_weight:
                            res += ring(wp_rc)
                        # NOT ACTUALLY UNIQUE RC
                        
    return res


if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult.abc import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing, RootTableau, Sx, CrystalGraphTensor, uncode, DSx
    # from schubmultutils.perm_utils import artin_sequences

    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)
    # dom_perms = {perm.minimal_dominant_above() for perm in perms}
    # THIS WORKS NOTE
    # for dom in dom_perms:
    #     for v in perms:
    #         prod = Sx(dom) * Sx(v)
    #         for w in prod:
    #             print(f"{w.trimcode}")
    #             coeff = 0
    #             for dom_rc in RCGraph.all_rc_graphs(dom, n-1):
    #                 for v_rc in RCGraph.all_rc_graphs(v, n-1):
    #                     dp_ret = v_rc.dualpieri(dom, w)
    #                     if len(dp_ret) > 0:
    #                         coeff += 1

    # TRY w0 only
    
    # magic_tensors = {}
    # for dom_perm in dom_perms:
    #     for v in perms:        
    #         magic_tensors[(v, dom_perm)] = set()
    #         prod = Sx(dom_perm) * Sx(v)
    #         for w in prod:
    #             print(f"{w.trimcode}")
    #             coeff = 0
    #             for dom_rc in RCGraph.all_rc_graphs(dom_perm, n):
    #                 for v_rc in RCGraph.all_rc_graphs(v, n):
    #                     dp_ret = v_rc.dualpieri(dom_perm, w)
    #                     if len(dp_ret) > 0:
    #                         magic_tensors[(v, dom_perm)].add(CrystalGraphTensor(dom_rc, v_rc))
    # BIJECTIVE RULE WHERE v can be made dominant within u's ascents
    from schubmult import RCGraphRing
    rc_ring = RCGraphRing()
    lr_rcs = {}
    for u in perms:
        for v in perms:
            
            if not dominifiable(u, v):
                continue

            print(f"Testing u={u.trimcode}, v={v.trimcode}")
            dom_perm = v.minimal_dominant_above()
            if v != dom_perm:
                print("TESTING DOMINANT ONLY TEMP")
                continue
            tensor_hw = set()
            # dom_perm = v
            prod0 = Sx(u) * Sx(dom_perm)
            prod = Sx(u) * Sx(v)
            diff_elem = (~v) * (dom_perm)
            sanity_prod = Sx.zero
            rc_prod = crystal_dom_product(RCGraph.principal_rc(dom_perm, n-1), u)
            pretty_print(rc_prod)
            sanity_prod = sum([coeff * Sx(rc.perm) for rc, coeff in rc_prod.items()])

            # test_sum = S.Zero
            # crystals = set()
            # all_tensor_hw = set()
            # u_rcs = RCGraph.all_rc_graphs(u, n-1)
            # for w in prod0:
            #     for dom_rc in RCGraph.all_rc_graphs(dom_perm, n-1):
            #         # construct_lr = True
            #         # if (u, dom_perm, w) in lr_rcs:
            #         #     u_rcs = lr_rcs[(u, dom_perm, w)]
            #         #     construct_lr = False
            #         # else:
                    
            #         #     lr_rcs[(u, dom_perm, w)] = set()
                        
            #         for u_rc in u_rcs:
            #             all_tensor_hw.add(CrystalGraphTensor(dom_rc, u_rc).to_highest_weight()[0])
            #             dp_ret = u_rc.dualpieri(dom_perm, w)
            #             if len(dp_ret) > 0:
            #                 # lr_rcs[(u, dom_perm, w)].add(u_rc)
            #                 pants =  CrystalGraphTensor(dom_rc, u_rc).to_highest_weight()[0]
            #                 if pants not in crystals:   
            #                     crystals.add(pants)
            #                     wt = pants.to_lowest_weight()[0].crystal_weight
            #                     wp_rcs = [rc for rc in RCGraph.all_rc_graphs(w, n-1, weight=wt) if rc.is_lowest_weight]
            #                     assert len(wp_rcs) == 1
            #                     wp_rc = wp_rcs[0]
            #                     try:
            #                         assert wt == wp_rc.length_vector
            #                     except AssertionError:
            #                         print("Weight mismatch detected!")
            #                         pretty_print(pants)
            #                         print(f"Expected weight: {RCGraph.principal_rc(w, n-1).length_vector}, got {wt}")
            #                         pretty_print(wp_rc)
            #                         assert False

            #                 else:
            #                     print("Duplicate crystal detected!")
            #                     pretty_print(pants)
            #                     assert False
            #                 down_w = w*(~diff_elem)
            #                 if down_w.inv == w.inv - diff_elem.inv:
            #                     sanity_prod += Sx(down_w)        
            #                 # print("For fun:")
            #                 # print(f"{u_rc.dualpieri(v, down_w)=}")
            try:
                assert (prod - sanity_prod).expand() == S.Zero
                print("Verified")
            except AssertionError:
                print("Mismatch detected!")
                print("Expected:")
                print(prod)
                print("Got:")
                print(sanity_prod)
                pretty_print(rc_prod)
                raise
            # print("Missed tensor crystals:")
            # for atw in all_tensor_hw:
            #     if atw not in crystals:
            #         pretty_print(atw)
                # for rc in RCGraph.all_hw_rcs(v, n-1):
                #     dom_rc = RCGraph.principal_rc(dom, n-1)
                #     tensor_hw = set()
                #     for v_rc in rc.full_crystal:
                #         tensor = CrystalGraphTensor(dom_rc, v_rc)
                #         tensor_hw.add(tensor.to_highest_weight()[0])
                    
                #     for thw in tensor_hw:
                #         for indiv_tensor in thw.full_crystal:
                #             dp_ret = indiv_tensor.factors[1].dualpieri(dom, w)
                #             if len(dp_ret) > 0:
                #                 # check if crystal isomorphic
                #                 w_rcs = RCGraph.all_rc_graphs(w, n-1, weight=indiv_tensor.crystal_weight)
                #                 assert len(w_rcs) > 0, f"No RC graphs for {w} with weight {indiv_tensor.crystal_weight}"
                #                 found_isomorphic = False
                #                 for w_rc in w_rcs:
                #                     if is_isomorphic_crystal(indiv_tensor, w_rc):
                #                         found_isomorphic = True
                #                         break
                #                 assert found_isomorphic, f"No isomorphic crystal found for tensor {indiv_tensor} and w {w}"
                #                 print("Successfully verified isomorphic crystal for tensor and w")
                #                 assert len(w_rcs) == 1
                #                 print("This crystal is unique")
                #                coeff += 1
                # print(f"Coefficient for {dom} * {v} to {w} is {coeff}, should be {prod[w]}")
                # assert coeff == prod[w], f"Coefficient mismatch for {dom} * {v} to {w}: got {coeff}, expected {prod[w]}"
                # print("Verified")
    
    # After the dominifiable test, add:
    # print("\n" + "="*60)
    # print("Testing tensor decomposition in dominifiable cases")
    # print("="*60)
    # test_decomposition(dominifiable_case=True)
    # test_decomposition(dominifiable_case=False)

    # attempt
    # lr_rcs = {}
    # zob = [0 for _ in range(20)]
    # for k in range(1, n):
    #     for v in perms:
    #         dom_perm = uncode([1] * k)
    #         print(f"Testing k={k}, v={v.trimcode}")
    #         ring = DSx([]).ring
    #         prod = Sx(v) * ring(dom_perm)
    #         prod0 = ring(dom_perm) * ring(v)
    #         sanity_prod = Sx.zero
    #         for w in prod0:
    #             for dom_rc in RCGraph.all_rc_graphs(dom_perm, n-1):  
    #                 for v_rc in RCGraph.all_rc_graphs(v, n-1):
    #                     dp_ret = v_rc.dualpieri(dom_perm, w)
    #                     if len(dp_ret) > 0:
    #                         #assert len(dp_ret) == 1
    #                         for pr in dp_ret:

    #                             # assert pr[1].perm.inv == 0
    #                             # pretty_print(pr)
    #                             sanity_prod += (ring.coeff_genset[1] ** (w.inv - dom_perm.inv - v.inv)) * Sx(w)
    #         try:
    #             assert (prod.expand() - sanity_prod.expand()).expand() == S.Zero
    #             print("Verified")
    #         except AssertionError:
    #             print("Mismatch detected!")
    #             print("Expected:")
    #             print(f"{prod=}")
    #             print(prod.expand())
    #             print("Got:")
    #             print(f"{sanity_prod=}")
    #             print(sanity_prod.expand())
    #             print(f"{(prod.expand() - sanity_prod.expand()).expand()=}")
    #             raise
    
    # Add at the end:
    # print("\n" + "="*60)
    # print("TESTING INTERLEAVING HYPOTHESIS")
    # print("="*60)
    
    # # Test simple concatenation first (known to work for pairs)
    # print("\n1. Testing simple concatenation (u_first):")
    # concatenation_works = test_interleaving_lr_rule(pattern='u_first')
    
    # # Test alternating interleaving (the general case)
    # print("\n2. Testing alternating interleaving:")
    # alternating_works = test_interleaving_lr_rule(pattern='alternating')
    
    # if concatenation_works:
    #     print("\n✓ Simple concatenation works (as expected)")
    
    # if alternating_works:
    #     print("\n🎉 MIRACLE: Alternating interleaving works!")
    #     print("This gives a fully general bijective LR rule!")
    # else:
    #     print("\nAlternating interleaving doesn't work directly")
    #     print("Need to understand the pattern of failures...")