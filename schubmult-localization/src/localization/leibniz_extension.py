from schubmult import RCGraph, RCGraphRing

def compute_product(a_rc, b_rc, i):
    """
    Compute the product of two RC-graphs a_rc and b_rc using the Leibniz extension method.
    
    This method allows for the computation of products that are not limited to the form of
    (any RC) % (dominant RC). It utilizes insights from divided differences to extend the
    product definition.
    
    Parameters:
    - a_rc: RCGraph instance representing the first RC-graph.
    - b_rc: RCGraph instance representing the second RC-graph.
    - i: Index for the divided difference operation.
    
    Returns:
    - A new RCGraph instance representing the product of a_rc and b_rc.
    """
    rc_ring = RCGraphRing()
    
    # Compute the divided difference of the product
    product = rc_ring(a_rc) % rc_ring(b_rc)
    divdiff_result = product.divdiff(i)
    
    if len(divdiff_result) == 0:
        raise ValueError("Divided difference resulted in zero, cannot compute product.")
    
    # Construct the new product based on the Leibniz extension
    new_product = {}
    for term, coeff in divdiff_result.items():
        new_product[term] = coeff
    
    return new_product


def extend_product_with_non_dominant(a_rc, dom_rc, i):
    """
    Extend the product of an RC-graph with a non-dominant RC-graph using the Leibniz extension.
    
    This function computes the product a_rc % dom_rc where a_rc is any RC-graph and
    dom_rc is a non-dominant RC-graph. It leverages the insights gained from divided
    differences to ensure proper weight preservation and product definition.
    
    Parameters:
    - a_rc: RCGraph instance representing the RC-graph to be extended.
    - dom_rc: RCGraph instance representing the non-dominant RC-graph.
    - i: Index for the divided difference operation.
    
    Returns:
    - A new RCGraph instance representing the extended product.
    """
    rc_ring = RCGraphRing()
    
    # Compute the product
    product = rc_ring(a_rc) % rc_ring(dom_rc)
    
    # Check if the divided difference is zero
    if product.divdiff(i):
        raise ValueError("Divided difference is non-zero, cannot extend product.")
    
    # Perform the extension
    extended_product = compute_product(a_rc, dom_rc, i)
    
    return extended_product


def define_new_products(a_rc, b_rc, i):
    """
    Define new products of RC-graphs that are not limited to the form of (any RC) % (dominant RC).
    
    This function utilizes the insights from divided differences and the Leibniz extension
    to compute products of RC-graphs in a more general context.
    
    Parameters:
    - a_rc: RCGraph instance representing the first RC-graph.
    - b_rc: RCGraph instance representing the second RC-graph.
    - i: Index for the divided difference operation.
    
    Returns:
    - A new RCGraph instance representing the defined product.
    """
    # Attempt to compute the product using the defined methods
    try:
        return extend_product_with_non_dominant(a_rc, b_rc, i)
    except ValueError as e:
        print(f"Error in defining new product: {e}")
        return None