from schubmult import RCGraph, RCGraphRing

class WeightTracker:
    def __init__(self):
        self.weights = {}

    def track_weight(self, rc_graph):
        """
        Track the weight of a given RC-graph.
        """
        weight = rc_graph.crystal_weight
        self.weights[rc_graph.perm.trimcode] = weight

    def get_weight(self, perm):
        """
        Retrieve the weight associated with a given permutation.
        """
        return self.weights.get(perm.trimcode, None)

    def clear_weights(self):
        """
        Clear all tracked weights.
        """
        self.weights.clear()

    def display_weights(self):
        """
        Display all tracked weights.
        """
        for perm, weight in self.weights.items():
            print(f"Permutation: {perm}, Weight: {weight}")

def compute_product_with_weights(rc_a, rc_b):
    """
    Compute the product of two RC-graphs and track the resulting weight.
    """
    rc_ring = RCGraphRing()
    product = rc_ring(rc_a) % rc_ring(rc_b)
    
    weight_tracker = WeightTracker()
    weight_tracker.track_weight(rc_a)
    weight_tracker.track_weight(rc_b)
    
    for rc in product:
        weight_tracker.track_weight(rc)

    return product, weight_tracker