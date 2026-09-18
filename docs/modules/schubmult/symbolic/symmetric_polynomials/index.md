<a id="schubmult.symbolic.symmetric_polynomials"></a>

# schubmult.symbolic.symmetric\_polynomials

Symbolic (factorial) elementary and complete symmetric polynomials as SymPy functions.

``E(p, k, *vars)``/``e`` and ``H(p, k, *vars)``/``h`` are the elementary and complete symmetric
polynomials of degree ``p`` in the first ``k`` generators, kept unevaluated so Schubert
expansions can be written in the SEM basis; the ``Factorial*`` variants carry a second
variable set. `functions` holds canonicalization and variable-splitting utilities, and
`qelem_sym` the quantum elementary symmetric polynomials.

