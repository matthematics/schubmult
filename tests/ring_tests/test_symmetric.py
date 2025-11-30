def test_elem_sym():
    from schubmult.abc import x, y, E
    from schubmult.symbolic import expand,S
    elem = E(2,5,[x[1],x[4],x[3],x[5],x[7]],[y[3],y[1],y[4],y[2]])
    assert expand(elem - elem.split_out_vars([x[3],x[5]]), func=True) == S.Zero
    assert elem.xreplace({y[2]: x[4]}) == E(2,4,[x[1],x[3],x[5],x[7]],[y[3],y[1],y[4]])


def test_complete_sym():
    from schubmult.abc import x, y, H
    from schubmult.symbolic import expand, S
    elem = H(5,4,x[1],x[4],x[3],x[5],y[3],y[1],y[4],y[2],y[7],y[9],y[11],y[15])
    assert expand(elem - elem.split_out_vars([x[3],x[1]]), func=True) == S.Zero
    assert elem.xreplace({y[2]: x[4]}) == H(5,3,[x[1],x[3],x[5]],[y[3],y[1],y[4],y[7],y[9],y[11],y[15]])