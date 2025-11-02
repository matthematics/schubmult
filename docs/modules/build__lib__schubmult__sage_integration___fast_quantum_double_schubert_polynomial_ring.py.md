# build/lib/schubmult/sage_integration/_fast_quantum_double_schubert_polynomial_ring.py



## FastQuantumDoubleSchubertPolynomialRing(R, num_vars, varname1, varname2, q_varname)

Return the FastQuantumDoubleSchubert polynomial ring over ``R`` on the X basis.

This is the basis made of the FastQuantumDoubleSchubert polynomials.

EXAMPLES::

        sage: X = FastQuantumDoubleSchubertPolynomialRing(ZZ); X
        Schubert polynomial ring with X basis over Integer Ring
        sage: TestSuite(X).run()
        sage: X(1)
        X[1]
        sage: X([1,2,3])*X([2,1,3])
        X[2, 1]
        sage: X([2,1,3])*X([2,1,3])
        X[3, 1, 2]
        sage: X([2,1,3])+X([3,1,2,4])
        X[2, 1] + X[3, 1, 2]
        sage: a = X([2,1,3])+X([3,1,2,4])
        sage: a^2
        X[3, 1, 2] + 2*X[4, 1, 2, 3] + X[5, 1, 2, 3, 4]

## class FastQuantumDoubleSchubertPolynomial_class(CombinatorialFreeModule.Element)



## class FastQuantumDoubleSchubertPolynomialRing_xbasis(CombinatorialFreeModule)



