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
            print(f"  Component is NOT Demazure (size {len(component)})")
            # Debug: how many highest/lowest weights?
            hw_count = sum(1 for e in component if all(
                e.raising_operator(i) is None 
                for i in range(1, e.crystal_length())
            ))
            lw_count = sum(1 for e in component if all(
                e.lowering_operator(i) is None 
                for i in range(1, e.crystal_length())
            ))
            print(f"    Highest weights: {hw_count}")
            print(f"    Lowest weights: {lw_count}")
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
                        print(f"\n!!! NON-DEMAZURE COMPONENT FOUND !!!")
                        print(f"u={u.trimcode}, v={v.trimcode}")
    
    print(f"\n=== Summary ===")
    print(f"Decomposes into Demazure: {decompose_count}")
    print(f"Has non-Demazure components: {dont_decompose_count}")
    
    if len(non_demazure_examples) > 0:
        print(f"\nFound {len(non_demazure_examples)} cases with non-Demazure components!")
    else:
        print("\nAll tensor products decompose into Demazure crystals!")
    
    return non_demazure_examples


if __name__ == "__main__":
    # test module functionality

    import itertools
    import sys

    from schubmult.abc import x
    from symengine import S
    from sympy import pretty_print

    from schubmult import Permutation, RCGraph, RCGraphRing, RootTableau, Sx, CrystalGraphTensor
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
    if False:
        lr_rcs = {}
        for u in perms:
            for v in perms:
                
                if not dominifiable(u, v):
                    continue
                print(f"Testing u={u.trimcode}, v={v.trimcode}")
                dom_perm = v.minimal_dominant_above()
                prod0 = Sx(u) * Sx(dom_perm)
                diff_elem = (~v) * (dom_perm)
                prod = Sx(u) * Sx(v)
                sanity_prod = Sx.zero
                for w in prod0:
                    for dom_rc in RCGraph.all_rc_graphs(dom_perm, n-1):
                        construct_lr = True
                        if (u, dom_perm, w) in lr_rcs:
                            u_rcs = lr_rcs[(u, dom_perm, w)]
                            construct_lr = False
                        else:
                            u_rcs = RCGraph.all_rc_graphs(u, n-1)
                            lr_rcs[(u, dom_perm, w)] = set()
                            
                        for u_rc in u_rcs:
                            if construct_lr:
                                dp_ret = u_rc.dualpieri(dom_perm, w)
                                if len(dp_ret) > 0:
                                    lr_rcs[(u, dom_perm, w)].add(u_rc)
                                else:
                                    continue
                            down_w = w*(~diff_elem)
                            if down_w.inv == w.inv - diff_elem.inv:
                                sanity_prod += Sx(down_w)
                try:
                    assert (prod - sanity_prod).expand() == S.Zero
                    print("Verified")
                except AssertionError:
                    print("Mismatch detected!")
                    print("Expected:")
                    print(prod)
                    print("Got:")
                    print(sanity_prod)
                    raise

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
    print("Testing tensor decomposition in dominifiable cases")
    print("="*60)
    #test_decomposition(dominifiable_case=True)
    test_decomposition(dominifiable_case=False)

