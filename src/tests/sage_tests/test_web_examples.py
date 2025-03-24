
def test_basic_sage_example():
    from sage.all import ZZ
    from schubmult.sage_integration import FastSchubertPolynomialRing, FastDoubleSchubertPolynomialRing, FastQuantumSchubertPolynomialRing
    SingleRing = FastSchubertPolynomialRing(ZZ, 100, "x")
    assert str(SingleRing([3,4,1,2])) == "Sx([3, 4, 1, 2])"
    assert SingleRing([3,4,1,2]) * SingleRing([5,1,4,2,3]) == SingleRing([7, 3, 4, 1, 2, 5, 6]) + SingleRing([7, 4, 2, 1, 3, 5, 6]) + SingleRing([7, 5, 1, 2, 3, 4, 6])
    DoubleRing = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
    assert str(DoubleRing([3,4,1,2], "z")) == "DSx([3, 4, 1, 2], 'z')"
    pr = DoubleRing._base_polynomial_ring
    assert DoubleRing([3,4,1,2])*DoubleRing([1,4,2,3]) == pr("y1*y2-y1*y4-y2*y4+y4^2")*DoubleRing([3, 4, 1, 2], 'y') + pr("-y1-y2+y4+y5")*DoubleRing([3, 5, 1, 2, 4], 'y') + DoubleRing([3, 6, 1, 2, 4, 5], 'y')
    SingleRingQ = FastQuantumSchubertPolynomialRing(ZZ, 100, "x")
    prq = SingleRingQ._polynomial_ring
    assert SingleRingQ([2,3,1,4]).expand() == prq("x1*x2 + q1")
    assert DoubleRing([1,4,2,3],"z") * SingleRing([3,4,1,2]) == pr("z1^2*z4^2")*DoubleRing([1, 4, 2, 3], 'z') + pr("z1^2*z4+z1^2*z5")*DoubleRing([1, 5, 2, 3, 4], 'z') + pr("z1^2")*DoubleRing([1, 6, 2, 3, 4, 5], 'z') + pr("z1*z4^2+z2*z4^2")*DoubleRing([2, 4, 1, 3], 'z') + pr("z1*z4+z2*z4+z1*z5+z2*z5")*DoubleRing([2, 5, 1, 3, 4], 'z') + pr("z1+z2")*DoubleRing([2, 6, 1, 3, 4, 5], 'z') + pr("z4^2")*DoubleRing([3, 4, 1, 2], 'z') + pr("z4+z5")*DoubleRing([3, 5, 1, 2, 4], 'z') + DoubleRing([3, 6, 1, 2, 4, 5], 'z')
    assert SingleRingQ([2,3,1,4]) * SingleRing([4,1,3,2]) == prq("-2*q1^2*q2+q1*q2*q3")*SingleRingQ([1]) + prq("q1*q2")*SingleRingQ([1, 3, 4, 2]) + prq("-q1^2")*SingleRingQ([2, 3, 1]) + prq("q1")*SingleRingQ([2, 4, 3, 1]) + prq("-q1*q2")*SingleRingQ([3, 1, 2]) + prq("q1")*SingleRingQ([3, 2, 4, 1]) + prq("q1")*SingleRingQ([3, 4, 1, 2]) + prq("-q1")*SingleRingQ([4, 2, 1, 3]) + SingleRingQ([5, 2, 3, 1, 4]) + SingleRingQ([5, 3, 1, 2, 4])
    assert SingleRing([1,3,2]) - pr("x1") - pr("x2") == 0
    DoubleRing.set_coproduct_indices((1,3))
    cprd = DoubleRing([4,1,5,2,3], "z").coproduct()
    assert str(cprd) == "(y1^2-y1*z2-y1*z3+z2*z3)*DSx([4, 1, 2, 3], 'z') # DSx([1], 'y') + (y1+y2-z2-z3)*DSx([4, 1, 2, 3], 'z') # DSx([2, 1], 'y') + DSx([4, 1, 2, 3], 'z') # DSx([3, 1, 2], 'y') + (y1-z3)*DSx([4, 2, 1, 3], 'z') # DSx([1], 'y') + DSx([4, 2, 1, 3], 'z') # DSx([2, 1], 'y') + DSx([4, 3, 1, 2], 'z') # DSx([1], 'y')"
