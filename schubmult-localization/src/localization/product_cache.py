from collections import defaultdict

class ProductCache:
    def __init__(self):
        self.cache = defaultdict(dict)

    def get_product(self, rc_a, rc_b):
        """
        Retrieve the cached product of two RC-graphs if it exists.
        """
        return self.cache[rc_a][rc_b] if rc_b in self.cache[rc_a] else None

    def set_product(self, rc_a, rc_b, product):
        """
        Cache the product of two RC-graphs.
        """
        self.cache[rc_a][rc_b] = product

    def clear_cache(self):
        """
        Clear the entire cache.
        """
        self.cache.clear()

    def compute_and_cache_product(self, rc_a, rc_b, product_function):
        """
        Compute the product of two RC-graphs using the provided product function,
        cache the result, and return it.
        """
        product = self.get_product(rc_a, rc_b)
        if product is None:
            product = product_function(rc_a, rc_b)
            self.set_product(rc_a, rc_b, product)
        return product

    def __repr__(self):
        return f"ProductCache({self.cache})"