from schubmult.localization.leibniz_extension import compute_leibniz_extension
from schubmult.products.extended_products import compute_extended_product
from schubmult.crystal.non_principal_actions import NonPrincipalRCGraph
from schubmult.products.product_cache import ProductCache
import unittest

class TestExtendedProducts(unittest.TestCase):

    def setUp(self):
        self.cache = ProductCache()

    def test_extended_product_with_dominant_rc(self):
        # Test case for computing extended product with a dominant RC
        dominant_rc = NonPrincipalRCGraph('dominant_permutation')
        non_dominant_rc = NonPrincipalRCGraph('non_dominant_permutation')
        
        result = compute_extended_product(non_dominant_rc, dominant_rc)
        
        # Assert expected properties of the result
        self.assertIsNotNone(result)
        self.assertTrue(result.has_valid_structure())

    def test_extended_product_with_non_dominant_rc(self):
        # Test case for computing extended product with a non-dominant RC
        rc1 = NonPrincipalRCGraph('non_dominant_permutation_1')
        rc2 = NonPrincipalRCGraph('non_dominant_permutation_2')
        
        result = compute_extended_product(rc1, rc2)
        
        # Assert expected properties of the result
        self.assertIsNotNone(result)
        self.assertTrue(result.has_valid_structure())

    def test_leibniz_extension_integration(self):
        # Test case to ensure Leibniz extension works with extended products
        rc = NonPrincipalRCGraph('some_permutation')
        dominant_rc = NonPrincipalRCGraph('dominant_permutation')
        
        leibniz_result = compute_leibniz_extension(rc, dominant_rc)
        extended_result = compute_extended_product(rc, dominant_rc)
        
        # Assert that the results are consistent
        self.assertEqual(leibniz_result, extended_result)

    def test_product_cache_usage(self):
        # Test case to verify caching mechanism
        rc1 = NonPrincipalRCGraph('permutation_1')
        rc2 = NonPrincipalRCGraph('permutation_2')
        
        product1 = self.cache.get_product(rc1, rc2)
        product2 = compute_extended_product(rc1, rc2)
        
        self.assertEqual(product1, product2)

if __name__ == '__main__':
    unittest.main()