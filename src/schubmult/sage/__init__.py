try:
    import sage  # noqa: F401
except ImportError:
    raise ImportError("SageMath is not installed. Please install SageMath to use this module.")

from ._fast_double_schubert_polynomial_ring import (
    FastDoubleSchubertPolynomial,
    FastDoubleSchubertPolynomialRing,
    FastDoubleSchubertPolynomialRing_base,
    FastQuantumDoubleSchubertPolynomialRing,
)
from ._fast_schubert_polynomial_ring import (
    FastQuantumSchubertPolynomialRing,
    FastSchubertPolynomial,
    FastSchubertPolynomialRing,
    FastSchubertPolynomialRing_base,
)

__all__ = [
    "FastDoubleSchubertPolynomial",
    "FastDoubleSchubertPolynomialRing",
    "FastDoubleSchubertPolynomialRing_base",
    "FastQuantumDoubleSchubertPolynomialRing",
    "FastQuantumSchubertPolynomialRing",
    "FastSchubertPolynomial",
    "FastSchubertPolynomialRing",
    "FastSchubertPolynomialRing_base",
]
