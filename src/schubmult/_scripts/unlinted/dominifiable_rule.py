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
    if len(tensor.full_crystal) != len(rc_graph.full_crystal):
        return False
    return True

# crystal level dominant product
# def crystal_dom_product(dom_rc, u_rc_hw):
#     from schubmult import RCGraphRing, Sx
#     u = u_rc_hw.perm
#     for w in cheat_prod:
#         for u_rc in u_rc_hw.full_crystal:
#             res += _rc_dom_product(dom_rc, u_rc, w, ring)
#     return res

def rc_dom_product(dom_rc, u_rc):
    ring = RCGraphRing()
    res = ring.zero
    cheat_prod = Sx(dom_rc.perm) * Sx(u_rc.perm)
    for w in cheat_prod:
        res += _rc_dom_product(dom_rc, u_rc, w, ring)
    return res

def rc_dom_crystal_product(dom_rc, u_rc_hw):
    # INSERTION WEIGHT TABLEAU
    rc_ring = RCGraphRing()
    tensor_hw_map = {}
    w_hw_map = {}
    n = len(dom_rc) + 1
    cheat_prod = Sx(dom_rc.perm) * Sx(u_rc_hw.perm)
    for u_rc_crystal in u_rc_hw.full_crystal:
        tensor_hw_map[u_rc_crystal] = CrystalGraphTensor(dom_rc, u_rc_crystal).to_highest_weight()[0]
        for w in cheat_prod:
            dp_ret = u_rc_crystal.dualpieri(dom_rc.perm, w)
            if len(dp_ret) > 0:
                pants =  tensor_hw_map[u_rc_crystal]
                wt = pants.to_lowest_weight()[0].crystal_weight
                wp_rcs = [rc for rc in RCGraph.all_rc_graphs(w, n-1, weight=wt) if rc.is_lowest_weight]
                wp_rc = wp_rcs[0]
                if wp_rc.to_highest_weight()[0].crystal_weight == pants.crystal_weight:
                    w_hw_map[pants] = wp_rc.to_highest_weight()[0]


def rc_dom_single_product(u_rc, v_rc):
    # INSERTION WEIGHT TABLEAU
    from schubmult import Plactic
    rc_ring = RCGraphRing()
    if v_rc.perm.minimal_dominant_above() == v_rc.perm:
        dom_rc = v_rc
        tensor_hw_map = {}
        w_hw_map = {}
        n = len(dom_rc) + 1
        cheat_prod = Sx(dom_rc.perm) * Sx(u_rc.perm)
        for u_rc_crystal in u_rc.full_crystal:
            tensor_hw_map[u_rc_crystal] = CrystalGraphTensor(dom_rc, u_rc_crystal).to_highest_weight()[0]
            for w in cheat_prod:
                dp_ret = u_rc_crystal.dualpieri(dom_rc.perm, w)
                if len(dp_ret) > 0:
                    pants =  tensor_hw_map[u_rc_crystal]
                    wt = pants.to_lowest_weight()[0].crystal_weight
                    wp_rcs = [rc for rc in RCGraph.all_rc_graphs(w, n-1, weight=wt) if rc.is_lowest_weight]
                    wp_rc = wp_rcs[0]
                    if wp_rc.to_highest_weight()[0].crystal_weight == pants.crystal_weight:
                        w_hw_map[pants] = wp_rc.to_highest_weight()[0]
        tensor =  tensor_hw_map[u_rc]
        tensor0 = CrystalGraphTensor(dom_rc, u_rc)
        tensor0_hw, raise_seq = tensor0.to_highest_weight()
        if tensor in w_hw_map:
            w_rc = w_hw_map[tensor]
            
            return rc_ring(w_rc.reverse_raise_seq(raise_seq))
        collected_rcs = {}
        for w in cheat_prod:
            for w_rc in RCGraph.all_hw_rcs(w, n-1):
                if w_rc not in w_hw_map.values():
                    collected_rcs[RootTableau.from_rc_graph(w_rc).weight_tableau] = w_rc

        tab1 = RootTableau.from_rc_graph(dom_rc).weight_tableau
        tab2 = RootTableau.from_rc_graph(tensor.factors[1]).weight_tableau

        total_tab = Plactic().rs_insert(*tab1.row_word, *tab2.row_word)
        try:
            return rc_ring(collected_rcs[total_tab].reverse_raise_seq(raise_seq))
        except KeyError:
            total_tab = Plactic().rs_insert(*tab2.row_word, *tab1.row_word)
            return rc_ring(collected_rcs[total_tab].reverse_raise_seq(raise_seq))
    # elif dominifiable(u_rc.perm, v_rc.perm):
    #     print("Warning: not fully working")
    #     actual_dom = v_rc.perm.minimal_dominant_above()
    #     diff_elem = (~v_rc.perm)*actual_dom
    #     dom_rc = RCGraph.principal_rc(actual_dom, len(v_rc))
    #     initial_prod = rc_dom_single_product(u_rc, dom_rc)
    #     candidates = initial_prod.divdiff_perm(diff_elem)
    #     try:
    #         return rc_ring(next(iter([rc for rc in candidates if rc.length_vector == CrystalGraphTensor(v_rc, u_rc).crystal_weight])))
    #     except StopIteration:
    #         print("No matching RC found in candidates")
    #         pretty_print(candidates)
    #         pretty_print(CrystalGraphTensor(v_rc, u_rc))
    #         raise
        # tensor_dom = CrystalGraphTensor(dom_rc, u_rc)
        # tensor_dom_hw, _ = tensor_dom.to_highest_weight()
        # tensor_hw, raise_seq = CrystalGraphTensor(v_rc, tensor_dom_hw.factors[1]).to_highest_weight()
        # collected_rcs = {}
        # for w_rc in candidates.keys():
        #     w_rc_hw = w_rc.to_highest_weight()[0]
        #     collected_rcs[RootTableau.from_rc_graph(w_rc_hw).weight_tableau] = w_rc_hw
        # tab1 = RootTableau.from_rc_graph(tensor_hw.factors[0]).weight_tableau
        # tab2 = RootTableau.from_rc_graph(tensor_hw.factors[1]).weight_tableau

        # total_tab = Plactic().rs_insert(*tab1.row_word, *tab2.row_word)
        # try:
        #     return rc_ring(collected_rcs[total_tab].reverse_raise_seq(raise_seq))
        # except KeyError:
        #     total_tab = Plactic().rs_insert(*tab2.row_word, *tab1.row_word)
        #     return rc_ring(collected_rcs[total_tab].reverse_raise_seq(raise_seq))
    raise NotImplementedError("General case not implemented yet")
        # WEIGHT TABLEAU INSERT
    
    

    #tab1 = dom_rc.weight_tableau
    # res = rc_ring.zero
    # for w_rc in w_hw_map.values():
    #     # THIS MAY FAIL: WHICH CRYSTALS MISSING? WEIGHT TABLEAU?
    #     # if tensor_hw not in w_hw_map:
    #     #     tab0 = RootTableau.from_rc_graph(tensor_hw.factors[0])
    #     #     tab1 = RootTableau.from_rc_graph(tensor_hw.factors[1])
    #     #     full_weight_tableau = Plactic().rs_insert(*tab0.weight_tableau, *tab1.weight_tableau)     
    #         # WILL HAVE SAME SHAPE
    #     res += rc_ring.from_dict(dict.fromkeys(w_rc.full_crystal,1))
    #     # for u_rc_hw in RCGraph.all_hw_rcs(u_rc.perm, n-1):
    #     #     # INSERT DOMINANT RC INTO u_rc_hw
    # while len(set(tensor_hw_map.values())) != len(w_hw_map):
    #     oldvals = set(tensor_hw_map.values())
    #     for tensor in oldvals:
    #         if tensor in w_hw_map:
    #             continue
    #         keyset = set(res.keys())
    #         perms_hit = set()
    #         for w_rc in keyset:
    #             w_rc_hw = w_rc.to_highest_weight()[0]
    #             for rc_hw in RCGraph.all_hw_rcs(w_rc.perm, n-1):
    #                 # if rc_hw in rcs_hit:
    #                 #     continue
                    
    #                 if is_isomorphic_crystal(tensor, rc_hw):
    #                     # missing crystal
    #                     coeff = res[w_rc]
    #                     res += rc_ring.from_dict(dict.fromkeys(rc_hw.full_crystal,coeff))
    #                     w_hw_map[tensor] = rc_hw
    #                     break
    #                 # rcs_hit.add(rc_hw)
    #                 # break
    # return res

#def lr_ed_decomp(dom_rc, perm):

def recursive_try_product(u_rc, v_rc, ring):
    if len(u_rc) != len(v_rc):
        return ring.zero
    if u_rc.inv == 0:
        return ring(v_rc)
    if v_rc.inv == 0:
        return ring(u_rc)
    if v_rc.perm.minimal_dominant_above() == v_rc.perm:
        return ring(u_rc) % ring(v_rc)
    
    cheat_prod = Sx(v_rc.perm) * Sx(u_rc.perm)
    if len(u_rc) == 1:
        result = ring.zero
        for rc in cheat_prod:
            return rc_ring(rc)
    
    # if mid_u < mid_v:
    #     mid_v = mid_u
    # elif mid_v < mid_u:
    #     mid_u = mid_v
    # if mid_u == len(u_rc):
    #     raise NotImplementedError("Cannot split further")
    left_u, right_u = u_rc.vertical_cut(len(u_rc) - 1)
    left_v, right_v = v_rc.vertical_cut(len(u_rc) - 1)
    left_prod = recursive_try_product(left_u, left_v, ring)
    right_prod = recursive_try_product(right_u, right_v, ring)

    combine_attempt = left_prod * right_prod
    result = ring.zero
    for rc in combine_attempt.keys():
        good = True
        if rc.vertical_cut(len(u_rc) - 1) == (next(iter(left_prod.keys())), next(iter(right_prod.keys()))) and rc.perm in cheat_prod:
            result += rc_ring(rc)
            break
    # total_prod = ring.zero
    # for left_rc in left_prod.keys():
    #     for right_rc in right_prod.keys():
    #         total_prod += (ring(left_rc) % ring(right_rc)) * (left_prod[left_rc] * right_prod[right_rc])
    # return total_prod
    # first cut
    print("Got through:")
    return result
    #return None

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

    for perm in perms:
        for rc in RCGraph.all_rc_graphs(perm, n-1):
            pretty_print(rc)
            print(f"{rc.params=}")

    input()

    dom_perms = {perm.minimal_dominant_above() for perm in perms}
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
    
    
    # MULTIPLY RC GRAPHS
    # BIJECTIVE RULE WHERE v can be made dominant within u's ascents
    from schubmult import RCGraphRing
    rc_ring = RCGraphRing()
   
    # Test module action
    for w in perms:
        for dom1, dom2 in itertools.product(dom_perms, dom_perms):
            dom1_rc = RCGraph.principal_rc(dom1, n-1)
            dom2_rc = RCGraph.principal_rc(dom2, n-1)

            for w_rc in RCGraph.all_rc_graphs(w, n-1):
                res1 =  rc_ring(w_rc) % (rc_ring(dom1_rc) % rc_ring(dom2_rc))
                res2 = (rc_ring(w_rc) % rc_ring(dom1_rc)) % rc_ring(dom2_rc)
                diff = res1 - res2
                assert all(coeff == 0 for coeff in diff.values())
            # cheat_prod = Sx(v) * Sx(u)
            # # CRYSTAL MONK
            # # dominifiable cauchy
            # # if dominifiable(u, v):
            # #     print(f"\n=== Testing dominifiable pair u={u.trimcode}, v={v.trimcode} ===")
            # # result = rc_ring.zero
            # # dom_rc = RCGraph.principal_rc(v.minimal_dominant_above(), n-1)
            # result = rc_ring.zero
            # for u_rc in RCGraph.all_rc_graphs(u, n-1):
            #     for v_rc in RCGraph.all_rc_graphs(v, n-1):
            #         # print("ATTEMPT")
            #         # pretty_print(CrystalGraphTensor(u_rc, v_rc))
            #         result += recursive_try_product(u_rc, v_rc, rc_ring)
            #         # except Exception as e:
            #         #     import traceback
            #         #     print("FAILED TO COMPUTE PRODUCT:")
            #         #     traceback.print_exc(file=sys.stdout)
            # #assert 
            #             #print(f"FAILED TO COMPUTE PRODUCT: {e}")
            # #         if v_rc.perm.minimal_dominant_above() != v_rc.perm:
            # #             continue
            # #         dom_rc = v_rc
            # #     result += (rc_ring(u_rc) % rc_ring(dom_rc)).divdiff_perm((~v)*dom_rc.perm)
            # #         # diff_elem = (~v)*v.minimal_dominant_above()
            # #         # final_prod = prod.divdiff_perm(diff_elem)
            # #         # cheat_prod = Sx(v) * Sx(u)
            # #         # total_count = sum(cheat_prod[w] for w in final_prod.keys())
            # #         # print(f"  u_rc perm={u_rc.perm.trimcode}, v_rc perm={v_rc.perm.trimcode} => count={total_count}")
            # # pretty_print(result)
            # test_prod = Sx(sum([coeff * rc.polyvalue(Sx.genset) for rc, coeff in result.items()]))
            # assert (test_prod - cheat_prod).expand() == S.Zero
            # print(f"Successfully computed product for u={u.trimcode}, v={v.trimcode}")
            # pretty_print(result)
            # print(f"Comparison: {cheat_prod}")


            # THEOREM: There is an action of the ring of dominant RC graphs on the module
            # of all RC graphs, and there is a surjective weight preserving homomorphism of coalgebras
            # and homomorphism of DRC-modules
            # from the module of all RC graphs to the polynomial ring that coincides with the action
            # and commutes with divided difference operators on the Schubert submodule.