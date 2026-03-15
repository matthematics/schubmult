from schubmult import *

if __name__ == "__main__":
    # from schubmult.rings.polynomial_algebra.anti_schubert_poly_basis import AntiSchubertPolyBasis
    # from schubmult.rings.polynomial_algebra import *

    # ASchub = PolynomialAlgebra(AntiSchubertPolyBasis())

    # perm = Permutation(uncode([0,5,2]))

    # elem = ASchub(perm, 3).change_basis(ASchub._basis.monomial_basis).change_basis(SchubertPolyBasis()).change_basis(KeyPolyBasis(Sx.genset))
    # assert Sx.from_expr(elem.expand()) == Sx.from_expr(ASchub(perm,3).expand())
    # #change_basis(SchubertPolyBasis())
    # print(elem)

    rc = RCGraph.random_rc_graph(uncode([0, 4, 2, 0, 5]))
    grass_rc = RCGraph.random_rc_graph(uncode([0, 2, 3, 3, 3]))
    grass_rc2 = RCGraph.random_rc_graph(uncode([0, 1, 2, 3, 4]))

    #print(rc.squash_product(grass_rc).squash_product(grass_rc2))

    anti_rc = AntiRCGraph.from_rc_graph(rc)
    anti_grass_rc = AntiRCGraph.from_rc_graph(grass_rc)
    anti_grass_rc2 = AntiRCGraph.from_rc_graph(grass_rc2)

    arc = rc.squash_product(grass_rc).squash_product(grass_rc2)
    print(arc)
    print(arc.length_vector)
    print(arc.perm.descents())
    prod = Sx(rc.perm) * Sx(grass_rc.perm) * Sx(grass_rc2.perm)
    assert prod.get(arc.perm, 0) != 0
    #print(rc.squash_product(grass_rc.squash_product(grass_rc2)))

    print("-----------------------")
    anti_arc = (anti_grass_rc.squash_product(anti_grass_rc2)).squash_product(anti_rc).to_rc_graph()
    print(anti_arc)
    print(anti_arc.length_vector)
    print(anti_arc.perm.descents())
    assert prod.get(anti_arc.perm, 0) != 0

    amidumb = rc.squash_product(grass_rc2.squash_product(grass_rc))
    assert amidumb == anti_arc
    print(amidumb)
    #print(anti_grass_rc.squash_product(anti_grass_rc2.squash_product(anti_rc)).to_rc_graph())
    
    
