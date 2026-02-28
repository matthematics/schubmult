from schubmult import *

if __name__ == "__main__":
    from schubmult.rings.free_algebra import *
    from schubmult.rings.polynomial_algebra import *
    from schubmult.utils.tuple_utils import *
    from schubmult.utils.perm_utils import weak_compositions
    ForestDual = FreeAlgebra(ForestBasis)
    from sympy import pretty_print
    import sys
    n = int(sys.argv[1])

    perms = Permutation.all_permutations(n)

     # We have two surjective ring homomorphisms from the RC graph ring to the ring of diagonal coinvariants: one sending an RC graph to its Schubert polynomial and one sending an RC graph to the monomial corresponding to its weight.
     # The fact that these are ring homomorphisms follows from the combinatorial definition of the product in the RC graph ring and the combinatorial definitions of the products in the Schubert basis and monomial basis of the diagonal coinvariant ring. The fact that they are surjective follows from the fact that every Schubert polynomial and every monomial can be realized as the image of some RC graph.
     # Moreover, these two maps are actually the same map, which can be seen from the fact that both maps send an RC graph to a homogeneous polynomial with leading term given by the monomial corresponding to its weight, and both maps are surjective onto a basis of the diagonal coinvariant ring.
     # Hence we have a single surjective ring homomorphism from the RC graph ring to the diagonal coinvariant ring that sends an RC graph to both its Schubert polynomial and its weight monomial, which are equal in the diagonal coinvariant ring.
     # By virtue of the fact that the diagonal coinvariant ring is free, there is a section of this map that sends a monomial in the diagonal coinvariant ring to a unique RC graph with that weight, which also has that Schubert polynomial as its image under this section. This gives us a bijective dictionary between RC graphs and monomials in the diagonal coinvariant ring (which are also Schubert polynomials) that respects the algebra structure.

    # Define the two maps from RC graphs to the diagonal coinvariant ring
    r = RCGraphRing()
    cauchy = (FA @ r).zero

    for cd in weak_compositions(n-1, n - 1):
        word_elem = ForestDual(*pad_tuple(cd, n - 1)).change_basis(WordBasis)
        forest_poly = ForestDual(*pad_tuple(cd, n - 1)).change_basis(SchubertBasis)
        for word, coeff0 in word_elem.items():
            for (perm, ln), coeff in forest_poly.items():
                cauchy += coeff * coeff0 * FA(*word) @ r.from_free_algebra_element(ASx(perm, ln))

    cauchy_tensor = (ForestDual @ r).zero

    for (word, rc), coeff in cauchy.items():
        cauchy_tensor += coeff * FA(*word).change_basis(ForestBasis) @ r(rc)

    print("Cauchy element in RC graph ring:")
    pretty_print(cauchy_tensor)