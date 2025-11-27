from schubmult.localization.leibniz_extension import compute_product
from schubmult.crystal.non_principal_actions import apply_non_principal_action

def check_product_consistency(rc_a, rc_b, product):
    """
    Verify the consistency of the product of two RC-graphs.
    
    Parameters:
    rc_a: RCGraph
        The first RC-graph.
    rc_b: RCGraph
        The second RC-graph.
    product: RCGraph
        The computed product of rc_a and rc_b.
    
    Returns:
    bool
        True if the product is consistent, False otherwise.
    """
    # Check if the product matches the expected properties
    expected_weight = tuple(map(sum, zip(rc_a.crystal_weight, rc_b.crystal_weight)))
    
    if product.crystal_weight != expected_weight:
        print(f"Weight mismatch: expected {expected_weight}, got {product.crystal_weight}")
        return False
    
    # Additional consistency checks can be added here
    # For example, checking the descents or structure of the product
    
    return True

def verify_extended_product(rc_a, rc_b):
    """
    Compute and verify the product of two RC-graphs, allowing for non-dominant cases.
    
    Parameters:
    rc_a: RCGraph
        The first RC-graph.
    rc_b: RCGraph
        The second RC-graph.
    
    Returns:
    bool
        True if the product is consistent, False otherwise.
    """
    # Compute the product using the Leibniz extension
    product = compute_product(rc_a, rc_b)
    
    # Check the consistency of the computed product
    return check_product_consistency(rc_a, rc_b, product)

def run_consistency_checks(rc_graphs):
    """
    Run consistency checks on a list of RC-graphs.
    
    Parameters:
    rc_graphs: list of RCGraph
        The list of RC-graphs to check.
    
    Returns:
    None
    """
    for i in range(len(rc_graphs)):
        for j in range(i, len(rc_graphs)):
            rc_a = rc_graphs[i]
            rc_b = rc_graphs[j]
            if not verify_extended_product(rc_a, rc_b):
                print(f"Inconsistent product for {rc_a.perm.trimcode} and {rc_b.perm.trimcode}")

# Example usage (to be removed or commented out in production)
# if __name__ == "__main__":
#     from schubmult import RCGraph
#     rc_graphs = [RCGraph.principal_rc(perm, n-1) for perm in Permutation.all_permutations(n)]
#     run_consistency_checks(rc_graphs)