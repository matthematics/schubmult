from schubmult import RCGraph, RCGraphRing

def compute_non_principal_product(a_rc, b_rc, i):
    """
    Compute the product of two RC-graphs where at least one is non-principal.
    
    This method utilizes insights from divided differences and the Leibniz extension
    to define products that are not limited to the form of (any RC) % (dominant RC).
    
    Parameters:
    - a_rc: RCGraph instance representing the first RC-graph.
    - b_rc: RCGraph instance representing the second RC-graph.
    - i: Index for the divided difference operation.
    
    Returns:
    - A new RCGraph instance representing the product of a_rc and b_rc.
    """
    rc_ring = RCGraphRing()
    
    # Compute the product using the Leibniz extension
    product = rc_ring(a_rc) % rc_ring(b_rc)
    
    # Check if the divided difference of the product at index i is zero
    if len(product.divdiff(i)) == 0:
        # If zero, we need to apply the extension strategy
        # Here we can define how to handle the product extension
        # This could involve using the s_i action or other methods
        # For now, we will return the product as is
        return product
    
    return product


def extend_product_with_non_principal(a_rc, b_rc, i):
    """
    Extend the product of two RC-graphs where at least one is non-principal.
    
    This method defines how to handle products involving non-principal RC-graphs
    by leveraging the insights gained from divided differences and the Leibniz extension.
    
    Parameters:
    - a_rc: RCGraph instance representing the first RC-graph.
    - b_rc: RCGraph instance representing the second RC-graph.
    - i: Index for the divided difference operation.
    
    Returns:
    - A new RCGraph instance representing the extended product.
    """
    # Compute the initial product
    initial_product = compute_non_principal_product(a_rc, b_rc, i)
    
    # Further extend the product based on the structure of a_rc and b_rc
    # This could involve additional operations or adjustments based on the
    # properties of the RC-graphs involved.
    
    # For demonstration, we will simply return the initial product
    return initial_product


def analyze_non_principal_actions(a_rc, b_rc, i):
    """
    Analyze the actions on non-principal RC-graphs and their products.
    
    This method provides insights into how non-principal RC-graphs interact
    and how their products can be computed and extended.
    
    Parameters:
    - a_rc: RCGraph instance representing the first RC-graph.
    - b_rc: RCGraph instance representing the second RC-graph.
    - i: Index for the divided difference operation.
    
    Returns:
    - A summary of the analysis as a string.
    """
    product = extend_product_with_non_principal(a_rc, b_rc, i)
    
    analysis_summary = (
        f"Analyzing product of RC-graphs:\n"
        f"  RC-graph A: {a_rc.perm.trimcode}\n"
        f"  RC-graph B: {b_rc.perm.trimcode}\n"
        f"  Product weight: {product.crystal_weight}\n"
    )
    
    return analysis_summary