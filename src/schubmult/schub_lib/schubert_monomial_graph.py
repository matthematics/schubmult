"""
Base class for Schubert monomial graph structures.

Provides a common interface for structures like RCGraph and BPD that represent
monomials in Schubert polynomials via grid-based diagrams.
"""

from abc import ABC, abstractmethod
from typing import Any

from schubmult.schub_lib.perm_lib import Permutation
from schubmult.symbolic import Expr


class SchubertMonomialGraph(ABC):
    """
    Abstract base class for Schubert monomial graph structures.

    This class provides a common interface for combinatorial objects that:
    - Represent monomials in Schubert polynomials
    - Have a grid/matrix-like structure with rows and columns
    - Correspond to a permutation
    - Have a weight (as a vector or tuple)

    Concrete implementations include:
    - RCGraph: Reduced-compatible graphs with crystal structure
    - BPD: Bumpless pipe dreams
    """

    @property
    @abstractmethod
    def perm(self) -> Permutation:
        """
        Return the permutation associated with this monomial graph.

        Returns:
            Permutation object
        """
        ...

    @property
    @abstractmethod
    def rows(self) -> int:
        """
        Number of rows in the grid representation.

        Returns:
            Number of rows
        """
        ...

    @property
    @abstractmethod
    def cols(self) -> int:
        """
        Number of columns in the grid representation.

        Returns:
            Number of columns
        """
        ...

    @property
    def width(self) -> int:
        """
        Width of the grid (alias for cols).

        Returns:
            Number of columns
        """
        return self.cols

    @property
    def height(self) -> int:
        """
        Height of the grid (alias for rows).

        Returns:
            Number of rows
        """
        return self.rows

    @property
    def permutation(self) -> Permutation:
        """
        Alias for perm property.

        Returns:
            Permutation object
        """
        return self.perm

    @abstractmethod
    def __getitem__(self, key) -> Any:
        """
        Access elements of the grid.

        Args:
            key: Index or tuple of indices

        Returns:
            Element(s) at the specified position
        """
        ...

    @abstractmethod
    def normalize(self) -> "SchubertMonomialGraph":
        """
        Return a normalized version of the monomial graph.

        Normalization may involve reordering or simplifying the structure
        while preserving its combinatorial properties.

        Returns:
            Normalized SchubertMonomialGraph object
        """
        ...

    @abstractmethod
    def polyvalue(self, x, y=None, **kwargs) -> Expr:
        """
        Compute the polynomial value represented by this monomial graph.

        Returns:
            Polynomial representation (e.g., as a SymPy expression)
        """
        ...

    @abstractmethod
    def shiftup(self, shift: int = 1) -> "SchubertMonomialGraph":
        """
        Shift up the monomial graph by a given amount.

        Args:
            shift: Amount to shift (default 1)

        Returns:
            Shifted monomial graph
        """
        ...

    @abstractmethod
    def right_zero_act(self) -> set["SchubertMonomialGraph"]:
        """
        Compute the right action of a zero (adding a row/column).

        Returns:
            Set of resulting monomial graphs
        """
        ...

    @abstractmethod
    def product(self, other: "SchubertMonomialGraph") -> dict["SchubertMonomialGraph", int]:
        """
        Compute the product of this monomial graph with another.

        This represents the Schubert polynomial multiplication at the monomial level.

        Args:
            other: Another monomial graph to multiply with

        Returns:
            Dictionary mapping result monomial graphs to their coefficients
        """
        ...
