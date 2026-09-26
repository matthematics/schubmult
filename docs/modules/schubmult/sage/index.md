<a id="schubmult.sage"></a>

# schubmult.sage

SageMath integration.

Importable only inside a Sage environment (``sage -python`` / a conda ``sagemath`` env with schubmult
installed). Provides Sage parents whose arithmetic is delegated to the schubmult kernels::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: X([3, 1, 2]) * X([2, 1])
    (y_2-y_0)*X_y[3, 1, 2] + X_y[4, 1, 2, 3]

