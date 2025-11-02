# src/schubmult/sage/_fast_double_schubert_polynomial_ring.py



## FastDoubleSchubertPolynomialRing(R, num_vars, base_variable_name, coeff_variable_names)

Wrapper function to return a double Schubert polynomial Ring

    Calls the _xbasis class to return a double or quantum double Schubert
    polynomial ring with the indicated base ring, number of variables,
    variable names (base variable, and then one or more sets of coefficient)
    variables, coproduct indices, code_display representation option, q-ring
    variable name, and whether the ring is quantum.

    Example call:

```python
X = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
X([2, 4, 3, 1]) + X([2, 1, 4, 3], "z")
```

Args:
        R (sage ring): The base ring
        num_vars (int): Cardinality of the sets of variables
        base_variable_name (str): Base variable name
        coeff_variable_names (str | tuple[str]): Coefficient variable name(s)
        indices (tuple[int], optional): Indicies of the variables to split on for the coproduct.
        code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.
        q_varname (str, optional): Variable name of the q-ring. Defaults to "q".
        is_quantum (bool, optional): Whether or not the ring is quantum. Defaults to False.

Returns:
        FastDoubleSchubertPolynomialRing_xbasis: Basis element generator of the ring

## FastQuantumDoubleSchubertPolynomialRing(R, num_vars, base_variable_name, coeff_variable_names)

Quantum double Schubert ring generator

Wraps FastDoubleSchubertPolynomialRing(), omitting indices and setting
is_quantum to True.

Args:
    R (sage ring): The base ring
    num_vars (int): Cardinality of the sets of variables
    base_variable_name (str): Base variable name
    coeff_variable_names (str | tuple[str]): Coefficient variable name(s)
    code_display (bool, optional): Whether to display the indices as the Lehmer code. Defaults to False.
    q_varname (str, optional): Variable name of the q-ring. Defaults to "q".

Returns:
    FastDoubleSchubertPolynomialRing_xbasis: Basis element generator of the quantum ring

## class FastDoubleSchubertPolynomial_class(CombinatorialFreeModule.Element)



## _double_schub_parser(passed)



## class FastDoubleSchubertPolynomialRing_xbasis(CombinatorialFreeModule)



