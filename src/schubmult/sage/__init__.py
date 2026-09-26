"""SageMath integration.

Importable only inside a Sage environment (``sage -python`` / a conda ``sagemath`` env with schubmult
installed). Provides Sage parents whose arithmetic is delegated to the schubmult kernels::

    sage: from schubmult.sage import DoubleSchubertPolynomialRing
    sage: X = DoubleSchubertPolynomialRing(QQ)
    sage: X([3, 1, 2]) * X([2, 1])
    (y_2-y_0)*X_y[3, 1, 2] + X_y[4, 1, 2, 3]

Rings: :func:`DoubleSchubertPolynomialRing`, :func:`QuantumSchubertPolynomialRing`,
:func:`QuantumDoubleSchubertPolynomialRing` (both with ``parabolic=``), :func:`GrothendieckPolynomialRing`
and :func:`DoubleGrothendieckPolynomialRing`.
"""

# Sage's submodules are not independently importable: entering the library at
# sage.categories.* first trips a circular import (``cannot import name Category``).
# ``sage.all`` is the supported bootstrap and is already loaded in any Sage session.
import sage.all  # isort: skip

from .double_schubert import DoubleSchubertPolynomialRing
from .grothendieck import DoubleGrothendieckPolynomialRing, GrothendieckPolynomialRing
from .quantum_schubert import QuantumDoubleSchubertPolynomialRing, QuantumSchubertPolynomialRing

__all__ = [
    "DoubleGrothendieckPolynomialRing",
    "DoubleSchubertPolynomialRing",
    "GrothendieckPolynomialRing",
    "QuantumDoubleSchubertPolynomialRing",
    "QuantumSchubertPolynomialRing",
]
