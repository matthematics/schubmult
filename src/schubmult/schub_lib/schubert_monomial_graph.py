"""
Base class for Schubert monomial graph structures.

Provides a common interface for structures like RCGraph and BPD that represent
monomials in Schubert polynomials via grid-based diagrams.
"""

from abc import ABC, abstractmethod

from schubmult.schub_lib.perm_lib import Permutation


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
        pass

    @property
    @abstractmethod
    def rows(self) -> int:
        """
        Number of rows in the grid representation.
        
        Returns:
            Number of rows
        """
        pass

    @property
    @abstractmethod
    def cols(self) -> int:
        """
        Number of columns in the grid representation.
        
        Returns:
            Number of columns
        """
        pass

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
    def __getitem__(self, key):
        """
        Access elements of the grid.
        
        Args:
            key: Index or tuple of indices
            
        Returns:
            Element(s) at the specified position
        """
        pass
