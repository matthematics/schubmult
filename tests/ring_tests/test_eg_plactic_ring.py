def test_to_from_rc_graph():
    from schubmult import EGPlacticRing, RCGraph, RCGraphRing
    r = RCGraphRing()
    e = EGPlacticRing()
    rc = RCGraph([(4,1),(3,2),(),(4,)])
    rc_elem = r(rc)
    assert rc_elem.almosteq(e.from_rc_graph(rc).to_rc_graph_ring_element())

def test_product_iso():
    from schubmult import EGPlacticRing, RCGraph, RCGraphRing
    r = RCGraphRing()
    e = EGPlacticRing()
    rc = RCGraph([(4,1),(3,2),(),(4,)])
    rc2 = RCGraph([(1,),(4,3,2),(4,),()])
    rc_elem = r(rc) * r(rc2)
    assert rc_elem.almosteq((e.from_rc_graph(rc) * e.from_rc_graph(rc2)).to_rc_graph_ring_element())