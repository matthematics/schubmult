# src/schubmult/sage/_fast_schubert_polynomial_ring.py



## FastSchubertPolynomialRing(R, num_vars, base_variable_name)

Wrapper function to return a double Schubert polynomial Ring

    Calls the _xbasis class to return a (quantum) Schubert
    polynomial ring with the indicated base ring, number of variables,
    variable name, coproduct indices, code_display representation option,
    q-ring variable name, and whether the ring is quantum.

    Example call:

```python
X = FastSchubertPolynomialRing(ZZ, 100, "x")
X([2, 4, 3, 1]) + X([2, 1, 4, 3])
```
This produces a sum of Schubert polynomials in the "x" variables. These will coerce
to any polynomial ring with variables with the same names as Schubert polynomials.

Args:
    R (Parent): The base ring
    num_vars (int): Cardinality of the sets of variables
    base_variable_name (str): Base variable name
    code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.
    q_varname (str, optional): Variable name of the q-ring. Defaults to "q".
    is_quantum (bool, optional): Whether or not the ring is quantum. Defaults to False.
    indices (tuple[int], optional): Indicies of the variables to split on for the coproduct.

Returns:
    FastSchubertPolynomialRing_xbasis: Element constructor of the ring

## FastQuantumSchubertPolynomialRing(R, num_vars, base_variable_name, q_varname, code_display)

Quantum Schubert ring generator

Wraps FastSchubertPolynomialRing(), omitting indices and setting
is_quantum to True.

Args:
    R (Parent): The base ring
    num_vars (int): Cardinality of the sets of variables
    base_variable_name (str): Base variable name
    q_varname (str, optional): Variable name of the q-ring. Defaults to "q".
    code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.


Returns:
    FastSchubertPolynomialRing_xbasis: Element constructor of the ring

## class FastSchubertPolynomial_class(CombinatorialFreeModule.Element)



## _single_schub_parser(passed)



## class FastSchubertPolynomialRing_xbasis(CombinatorialFreeModule)



