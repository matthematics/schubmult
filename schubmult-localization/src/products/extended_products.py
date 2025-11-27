from schubmult.localization.leibniz_extension import compute_leibniz_extension
from schubmult.crystal.non_principal_actions import apply_non_principal_action
from schubmult import RCGraph, RCGraphRing

def extended_product(a_rc, b_rc):
    """
    Compute the product of two RC-graphs, allowing for combinations
    that are not limited to the form (any RC) % (dominant RC).
    
    Parameters:
    - a_rc: RCGraph instance representing the first RC-graph.
    - b_rc: RCGraph instance representing the second RC-graph.
    
    Returns:
    - A new RCGraph instance representing the product of a_rc and b_rc.
    """
    rc_ring = RCGraphRing()
    
    # Step 1: Compute the Leibniz extension for the product
    leibniz_result = compute_leibniz_extension(a_rc, b_rc)
    
    # Step 2: Apply non-principal actions if necessary
    if not b_rc.perm.is_dominant:
        non_principal_result = apply_non_principal_action(leibniz_result, b_rc)
        return non_principal_result
    
    return leibniz_result

def compute_extended_products(rc_graphs):
    """
    Compute extended products for a list of RC-graphs.
    
    Parameters:
    - rc_graphs: List of RCGraph instances to compute products for.
    
    Returns:
    - A dictionary mapping pairs of RC-graphs to their computed products.
    """
    products = {}
    
    for i in range(len(rc_graphs)):
        for j in range(len(rc_graphs)):
            a_rc = rc_graphs[i]
            b_rc = rc_graphs[j]
            product = extended_product(a_rc, b_rc)
            products[(a_rc.perm.trimcode, b_rc.perm.trimcode)] = product
    
    return products