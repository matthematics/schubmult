<a id="schubmult.rings.free_algebra"></a>

# schubmult.rings.free\_algebra

Free algebra module providing multiple bases for Schubert calculus.

The core classes are :class:`FreeAlgebra` (the ring) and
:class:`FreeAlgebraElement` (its elements).  Elements can be expressed in
any of the available bases and converted between them via ``change_basis``.

Pre-built instances:
    - ``FA``: FreeAlgebra with WordBasis (default)
    - ``ASx``: FreeAlgebra with SchubertBasis
    - ``AGx``: FreeAlgebra with GrothendieckBasis
    - ``ADSx``: FreeAlgebra with double Schubert basis

Available bases:
    WordBasis, SchubertBasis, CompositionSchubertBasis, ElementaryBasis,
    ForestBasis, FundamentalSlideBasis, JBasis, JTBasis, KeyBasis,
    MonomialSlideBasis, NElementaryBasis, SchubertSchurBasis,
    SchurElementaryBasis, SeparatedDescentsBasis, ZBasis.

