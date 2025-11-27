from schubmult import RCGraph, RCGraphRing

def analyze_string_structure(rc_graph):
    """
    Analyze the structure of strings in the context of the given RC-graph.
    
    Parameters:
    rc_graph (RCGraph): The RC-graph to analyze.
    
    Returns:
    dict: A dictionary containing insights about the string structure.
    """
    insights = {
        'weight': rc_graph.crystal_weight,
        'descents': rc_graph.perm.descents(),
        'highest_weight': rc_graph.highest_weight(),
        'string_length': len(rc_graph.string()),
    }
    return insights


def compute_new_product(rc_a, rc_b):
    """
    Compute a new product of two RC-graphs that is not limited to the form of (any RC) % (dominant RC).
    
    Parameters:
    rc_a (RCGraph): The first RC-graph.
    rc_b (RCGraph): The second RC-graph.
    
    Returns:
    dict: A dictionary representing the new product of the two RC-graphs.
    """
    rc_ring = RCGraphRing()
    
    # Compute the product using the ring structure
    product = rc_ring(rc_a) % rc_ring(rc_b)
    
    # Analyze the resulting product
    product_insights = analyze_string_structure(product)
    
    return {
        'product': product,
        'insights': product_insights,
    }


def explore_non_principal_products(rc_graphs):
    """
    Explore products involving non-principal RC-graphs.
    
    Parameters:
    rc_graphs (list): A list of RC-graphs to explore.
    
    Returns:
    list: A list of computed products and their insights.
    """
    results = []
    
    for i, rc_a in enumerate(rc_graphs):
        for rc_b in rc_graphs[i + 1:]:
            result = compute_new_product(rc_a, rc_b)
            results.append(result)
    
    return results