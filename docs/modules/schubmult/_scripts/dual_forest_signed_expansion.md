<a id="schubmult._scripts.dual_forest_signed_expansion"></a>

# schubmult.\_scripts.dual\_forest\_signed\_expansion

CLI: dual forest polynomial computations from arXiv:2306.10939 (Nadeau-Tewari).

Two formulas are supported, both expressed as elements of
``FreeAlgebra(WordBasis)`` keyed by composition (exponent) tuples:

  * ``--mode dual`` (default):  signed lower-ideal expansion, equal to P_F by
    Theorem 4.1::

        P_F = sum_{lower ideals L of F} (-1)^|L|
                  * sum_{L-compatible kappa} x^kappa

  * ``--mode tilde``:  the unsigned dual forest polynomial of eq. (4.1)::

        ~P_F = sum_{internal(F)-compatible kappa} x^kappa

  * ``--mode both``:  print both.

