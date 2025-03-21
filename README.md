# schubmult

## Program and package for rapid computation of Littlewood-Richardson coefficients of Schubert polynomials, with optional Sage integration

The main purpose of this python package is for executing scripts to compute coefficients of products of various types of Schubert polynomials. Coproducts can also be computed, as well as substitution of commuting difference operators for quantum double Schubert polynomials. Quantum multiplication also has parabolic subgroup support, computed via the Peterson-Woodward comparison theorem. **Note that except for quantum Schubert polynomial multiplication with the --basic-pieri option, the methodology for quantum/quantum double Schubert polynomials is conjectural at this time.**


## Basic script command lines, one-line notation

```bash
schubmult_py 1 2 4 9 11 6 8 12 3 5 7 10 - 6 8 1 2 3 4 7 10 12 14 5 9 11 13  
schubmult_double 1 3 4 6 2 5 - 2 1 5 7 3 4 6  
schubmult_yz 1 3 4 6 2 5 - 2 1 5 7 3 4 6 --display-positive
schubmult_q 5 1 4 3 2 - 5 1 3 4 2
schubmult_q_double 5 1 4 3 2 - 5 1 3 4 2
schubmult_q_yz 5 1 4 3 2 - 2 5 1 3 4 --display-positive
```

## Using the Lehmer code

The same execution with the Lehmer code:

```bash
schubmult_py --code 0 0 1 5 6 2 3 4 - 5 6 0 0 0 0 1 2 3 4
schubmult_double --code 0 1 1 2 - 1 0 2 3
schubmult_yz --code 0 1 1 2 - 1 0 2 3 --display-positive
schubmult_q --code 4 0 2 1 - 4 0 1 1
schubmult_q_double --code 4 0 2 1 - 4 0 1 1
schubmult_q_yz --code 4 0 2 1 - 1 3 --display-positive
```

## Coproducts

Coproducts partition the variables of the polynomial ring and express the single polynomial as a sum of products of Schubert/double Schubert polynomials in the partitioned variables.

```bash
schubmult_py --coprod 1 3 5 7 2 4 6 - 2 4
schubmult_double --coprod 1 3 5 7 2 4 6 - 2 4
schubmult_yz --coprod 1 3 5 7 2 4 6 - 2 4 --display-positive
```
Equivalently with the Lehmer code:
```bash
schubmult_py --code --coprod 0 1 2 3 - 2 4
schubmult_double --code --coprod 0 1 2 3 - 2 4
schubmult_yz --code --coprod 0 1 2 3 - 2 4 --display-positive
```

## Quantum commuting difference operators

schubmult_q_yz has a feature for displaying the coefficients of the divided difference operators in the evaluation of the quantum double Schubert polynomials on the commuting difference operators of Fomin, Gelfand, and Postnikov. It is necessary to cap the value of n in the group S_n we are working in because as n increases the expression does not stabilize.
```bash
schubmult_q_yz --nil-hecke 6 --code 2 2 --display-positive
```

## Diplaying the result positively

The command line argument `--display-positive `is available in schubmult_yz and schubmult_q_yz, which displays the result positively (if possible, this is still only always possible conjecturally). It will fail and print out the offending case if it finds a counterexample. This is highly processor intensive.


Runtime will vary tremendously by case. The general problem is #P-hard. Though the result is always nonnegative (which at least is known for schubmult_py, schubmult_q, schubmult_double, and schubmult_q_double) and the problem is in GapP, it is not known to be in #P at this time.

schubmult_py is for multiplying ordinary Schubert polynomials. schubmult_yz is for multiplying double Schubert polynomials in different sets of coefficient variables (labeled y and z), and schubmult_double is for multiplying double Schubert polynomials in the same set of coefficient variables. Similarly, schubmult_q is for multiplying quantum Schubert polynomials, schubmult_q_double is for multiplying quantum double Schubert polynomials in the same set of coefficient variables, and schubmult_q_yz is for multiplying quantum double Schubert polynomials in different sets of coefficient variables, or in other words it computes the Gromov-Witten invariants, equivariant Gromov-Witten invariants, and (mixed?) equivariant Gromov-Witten invariants of the complete flag variety. All have the same command line syntax as schubmult, except when using the --code option. schubmult_double/schubmult_q_double display the result with nonnegative coefficients in terms of the negative simple roots (and the q variables), and schubmult_yz and schubmult_q_yz optionally display the result positively in terms of y_i-z_j (and q) with the --display-positive option.

schubmult_xx --coprod allows you to split (double) Schubert polynomials along certain indices (not available for quantum). It takes one permutation as an argument, followed by a dash -, then the set of indices you would like to split on. These coefficients are always nonnegative since they occur as product coefficients (this is actually how they are computed).

When imported as a python package, the relevant packages are schubmult.perm_lib, which has various permutation manipulation functions, and three modules that have functions of the same name (function name is "schubmult"): schubmult.schubmult_py, schubmult.schubmult_yz, schubmult.schubmult_double. Function takes a permutation dictionary (keys are tuples of ints, which must be trimmed permutations, and values are either integers or symengine values, which can also be integers) as well as a permutation as its second argument, which is the (double) Schubert polynomial to multiply by. Returns a dictionary of the same form with the coefficients.

```python
from schubmult.schubmult_yz import schubmult  
  
coeff_dict = schubmult({(1,3,4,6,2,5): 1},(2,1,5,7,3,4,6))  # outputs dictionary with results  
```


# Sage integration (as of version 1.5.0)

[SageMath](https://www.sagemath.org/) is a computer algebra system that, while wonderful, is monstrously large and only works on posix-based operating systems (including WSL VMs, so it is still usable on Windows). This is why Sage support is provided optionally in schubmult. The syntax to install the Sage dependencies is

```
pip install schubmult[sage]
```

This will install the [sagemath-standard](https://pypi.org/project/sagemath-standard/) python package in addition to the other dependencies. **Again, this only works on Linux, MacOS, or WSL.** To use with a currently installed SageMath distribution, use sage's python interpreter to install the package (the `[sage]` piece is not required in that case).

## Basic sage example

```python
sage: from schubmult.sage_integration import FastSchubertPolynomial, FastDoubleSchubertPolynomialRing, FastQuantumSchubertPolynomial, FastQuantumDoubleSchubertPolynomialRing
sage: SingleRing = FastSchubertPolynomialRing(ZZ, 100, "x")
sage: SingleRing([3,4,1,2])
Sx[3, 4, 1, 2]
sage: SingleRing([3,4,1,2]) * SingleRing([5,1,4,2,3])
Sx[7, 3, 4, 1, 2, 5, 6] + Sx[7, 4, 2, 1, 3, 5, 6] + Sx[7, 5, 1, 2, 3, 4, 6]
```

## Mixed variable (Molev-Sagan) type products

```python
sage: DoubleRing = FastDoubleSchubertPolynomialRing(ZZ, 100, "x", ("y", "z"))
sage: DoubleRing([3,4,1,2])*DoubleRing([1,4,2,3])
(y1*y2-y1*y4-y2*y4+y4^2)*Sx([3, 4, 1, 2], 'y') + (-y1-y2+y4+y5)*Sx([3, 5, 1, 2, 4], 'y') + Sx([3, 6, 1, 2, 4, 5], 'y')
sage: DoubleRing([3,4,1,2]) * DoubleRing([1,4,2,3],"z")
(y3^2+y3*y4+y4^2-y3*z1-y4*z1-y3*z2-y4*z2+z1*z2-y3*z3-y4*z3+z1*z3+z2*z3)*Sx([3, 4, 1, 2], 'y') + (y3+y4+y5-z1-z2-z3)*Sx([3, 5, 1, 2, 4], 'y') + Sx([3, 6, 1, 2, 4, 5], 'y')
```

## expand()

```python
sage: SingleRingQ = FastQuantumSchubertPolynomialRing(ZZ, 100, "x")
sage: SingleRingQ([2,3,1,4]).expand()
x1*x2 + q_1
```

## Coercion

Coercion was implemented as widely as possible.
```python
sage: DoubleRing([1,4,2,3],"z") * SingleRing([3,4,1,2])
z1^2*z4^2*Sx([1, 4, 2, 3], 'z') + (z1^2*z4+z1^2*z5)*Sx([1, 5, 2, 3, 4], 'z') + z1^2*Sx([1, 6, 2, 3, 4, 5], 'z') + (z1*z4^2+z2*z4^2)*Sx([2, 4, 1, 3], 'z') + (z1*z4+z2*z4+z1*z5+z2*z5)*Sx([2, 5, 1, 3, 4], 'z') + (z1+z2)*Sx([2, 6, 1, 3, 4, 5], 'z') + z4^2*Sx([3, 4, 1, 2], 'z') + (z4+z5)*Sx([3, 5, 1, 2, 4], 'z') + Sx([3, 6, 1, 2, 4, 5], 'z')
sage: SingleRingQ([2,3,1,4]) * SingleRing([4,1,3,2])
(-2*q_1^2*q_2+q_1*q_2*q_3)*QSx[1] + q_1*q_2*QSx[1, 3, 4, 2] + (-q_1^2)*QSx[2, 3, 1] + q_1*QSx[2, 4, 3, 1] + (-q_1*q_2)*QSx[3, 1, 2] + q_1*QSx[3, 2, 4, 1] + q_1*QSx[3, 4, 1, 2] + (-q_1)*QSx[4, 2, 1, 3] + QSx[5, 2, 3, 1, 4] + QSx[5, 3, 1, 2, 4]
sage: R.<x1, x2> = PolynomialRing(ZZ, 2)
sage: SingleRing([1,3,2]) - x1 - x2 == 0
True
```

## Coproducts

FastSchubertPolynomialRing and FastDoubleSchubertPolynomialRings are bialgebras and each element implements the `coproduct()` member function. `set_coproduct_indices()` on the base ring will determine the variables to partition on.
```ada
sage: DoubleRing.set_coproduct_indices((1,3))
sage: DoubleRing([4,1,5,2,3], "z").coproduct()
(y1^2-y1*z2-y1*z3+z2*z3)*Sx([4, 1, 2, 3], 'z') # Sx([1], 'y') + (y1+y2-z2-z3)*Sx([4, 1, 2, 3], 'z') # Sx([2, 1], 'y') + Sx([4, 1, 2, 3], 'z') # Sx([3, 1, 2], 'y') + (y1-z3)*Sx([4, 2, 1, 3], 'z') # Sx([1], 'y') + Sx([4, 2, 1, 3], 'z') # Sx([2, 1], 'y') + Sx([4, 3, 1, 2], 'z') # Sx([1], 'y')
```

## Demonstration of quantum double mixed products

```python
sage: QuantumDoubleRing([4,1,3,2])*QuantumDoubleRing([5,1,3,2,4], "z")
(q_1*q_2*q_3*y1^3+q_1*q_2*q_3*y1^2*y4+q_1*q_2*q_3*y1*y4^2+q_1*q_2*q_3*y4^3-q_1*q_2*q_3*y1^2*z1-q_1*q_2*q_3*y1*y4*z1-q_1*q_2*q_3*y4^2*z1-q_1*q_2*q_3*y1^2*z2-q_1*q_2*q_3*y1*y4*z2-q_1*q_2*q_3*y4^2*z2+q_1*q_2*q_3*y1*z1*z2+q_1*q_2*q_3*y4*z1*z2-q_1*q_2*q_3*y1^2*z3-q_1*q_2*q_3*y1*y4*z3-q_1*q_2*q_3*y4^2*z3+q_1*q_2*q_3*y1*z1*z3+q_1*q_2*q_3*y4*z1*z3+q_1*q_2*q_3*y1*z2*z3+q_1*q_2*q_3*y4*z2*z3-q_1*q_2*q_3*z1*z2*z3-q_1*q_2*q_3*y1^2*z4-q_1*q_2*q_3*y1*y4*z4-q_1*q_2*q_3*y4^2*z4+q_1*q_2*q_3*y1*z1*z4+q_1*q_2*q_3*y4*z1*z4+q_1*q_2*q_3*y1*z2*z4+q_1*q_2*q_3*y4*z2*z4-q_1*q_2*q_3*z1*z2*z4+q_1*q_2*q_3*y1*z3*z4+q_1*q_2*q_3*y4*z3*z4-q_1*q_2*q_3*z1*z3*z4-q_1*q_2*q_3*z2*z3*z4)*QSx([1], 'y') + (q_1*q_2*q_3*y1^2+q_1*q_2*q_3*y1*y4+q_1*q_2*q_3*y4^2+q_1*q_2*q_3*y1*y5+q_1*q_2*q_3*y4*y5+q_1*q_2*q_3*y5^2-q_1*q_2*q_3*y1*z1-q_1*q_2*q_3*y4*z1-q_1*q_2*q_3*y5*z1-q_1*q_2*q_3*y1*z2-q_1*q_2*q_3*y4*z2-q_1*q_2*q_3*y5*z2+q_1*q_2*q_3*z1*z2-q_1*q_2*q_3*y1*z3-q_1*q_2*q_3*y4*z3-q_1*q_2*q_3*y5*z3+q_1*q_2*q_3*z1*z3+q_1*q_2*q_3*z2*z3-q_1*q_2*q_3*y1*z4-q_1*q_2*q_3*y4*z4-q_1*q_2*q_3*y5*z4+q_1*q_2*q_3*z1*z4+q_1*q_2*q_3*z2*z4+q_1*q_2*q_3*z3*z4)*QSx([1, 2, 3, 5, 4], 'y') + (q_1*q_2*q_3*y1+q_1*q_2*q_3*y4+q_1*q_2*q_3*y5+q_1*q_2*q_3*y6-q_1*q_2*q_3*z1-q_1*q_2*q_3*z2-q_1*q_2*q_3*z3-q_1*q_2*q_3*z4)*QSx([1, 2, 3, 6, 4, 5], 'y') + q_1*q_2*q_3*QSx([1, 2, 3, 7, 4, 5, 6], 'y') + (q_1*q_3*y1^3+q_1*q_3*y1^2*y4+q_1*q_3*y1*y4^2+q_1*q_3*y4^3-q_1*q_3*y1^2*z1-q_1*q_3*y1*y4*z1-q_1*q_3*y4^2*z1-q_1*q_3*y1^2*z2-q_1*q_3*y1*y4*z2-q_1*q_3*y4^2*z2+q_1*q_3*y1*z1*z2+q_1*q_3*y4*z1*z2-q_1*q_3*y1^2*z3-q_1*q_3*y1*y4*z3-q_1*q_3*y4^2*z3+q_1*q_3*y1*z1*z3+q_1*q_3*y4*z1*z3+q_1*q_3*y1*z2*z3+q_1*q_3*y4*z2*z3-q_1*q_3*z1*z2*z3-q_1*q_3*y1^2*z4-q_1*q_3*y1*y4*z4-q_1*q_3*y4^2*z4+q_1*q_3*y1*z1*z4+q_1*q_3*y4*z1*z4+q_1*q_3*y1*z2*z4+q_1*q_3*y4*z2*z4-q_1*q_3*z1*z2*z4+q_1*q_3*y1*z3*z4+q_1*q_3*y4*z3*z4-q_1*q_3*z1*z3*z4-q_1*q_3*z2*z3*z4)*QSx([1, 4, 2, 3], 'y') + (q_1*y1^3*y3+q_1*y1^3*y4+q_1*y1^2*y3*y4+q_1*y1^2*y4^2+q_1*y1*y3*y4^2+q_1*y1*y4^3+q_1*y3*y4^3-q_1*y1^3*z1-q_1*y1^2*y3*z1-2*q_1*y1^2*y4*z1-q_1*y1*y3*y4*z1-2*q_1*y1*y4^2*z1-q_1*y3*y4^2*z1-q_1*y4^3*z1+q_1*y1^2*z1^2+q_1*y1*y4*z1^2+q_1*y4^2*z1^2-q_1*y1^3*z2-q_1*y1^2*y3*z2-2*q_1*y1^2*y4*z2-q_1*y1*y3*y4*z2-2*q_1*y1*y4^2*z2-q_1*y3*y4^2*z2-q_1*y4^3*z2+2*q_1*y1^2*z1*z2+q_1*y1*y3*z1*z2+3*q_1*y1*y4*z1*z2+q_1*y3*y4*z1*z2+2*q_1*y4^2*z1*z2-q_1*y1*z1^2*z2-q_1*y4*z1^2*z2+q_1*y1^2*z2^2+q_1*y1*y4*z2^2+q_1*y4^2*z2^2-q_1*y1*z1*z2^2-q_1*y4*z1*z2^2-q_1*y1^2*y3*z3-q_1*y1^2*y4*z3-q_1*y1*y3*y4*z3-q_1*y1*y4^2*z3-q_1*y3*y4^2*z3+q_1*y1^2*z1*z3+q_1*y1*y3*z1*z3+2*q_1*y1*y4*z1*z3+q_1*y3*y4*z1*z3+q_1*y4^2*z1*z3-q_1*y1*z1^2*z3-q_1*y4*z1^2*z3+q_1*y1^2*z2*z3+q_1*y1*y3*z2*z3+2*q_1*y1*y4*z2*z3+q_1*y3*y4*z2*z3+q_1*y4^2*z2*z3-2*q_1*y1*z1*z2*z3-q_1*y3*z1*z2*z3-2*q_1*y4*z1*z2*z3+q_1*z1^2*z2*z3-q_1*y1*z2^2*z3-q_1*y4*z2^2*z3+q_1*z1*z2^2*z3-q_1*y1^2*y3*z4-q_1*y1^2*y4*z4-q_1*y1*y3*y4*z4-q_1*y1*y4^2*z4-q_1*y3*y4^2*z4+q_1*y1^2*z1*z4+q_1*y1*y3*z1*z4+2*q_1*y1*y4*z1*z4+q_1*y3*y4*z1*z4+q_1*y4^2*z1*z4-q_1*y1*z1^2*z4-q_1*y4*z1^2*z4+q_1*y1^2*z2*z4+q_1*y1*y3*z2*z4+2*q_1*y1*y4*z2*z4+q_1*y3*y4*z2*z4+q_1*y4^2*z2*z4-2*q_1*y1*z1*z2*z4-q_1*y3*z1*z2*z4-2*q_1*y4*z1*z2*z4+q_1*z1^2*z2*z4-q_1*y1*z2^2*z4-q_1*y4*z2^2*z4+q_1*z1*z2^2*z4+q_1*y1*y3*z3*z4+q_1*y1*y4*z3*z4+q_1*y3*y4*z3*z4-q_1*y1*z1*z3*z4-q_1*y3*z1*z3*z4-q_1*y4*z1*z3*z4+q_1*z1^2*z3*z4-q_1*y1*z2*z3*z4-q_1*y3*z2*z3*z4-q_1*y4*z2*z3*z4+q_1*z1*z2*z3*z4+q_1*z2^2*z3*z4)*QSx([1, 4, 3, 2], 'y') + (q_1*y1^3+q_1*y1^2*y4+q_1*y1*y4^2+q_1*y4^3-q_1*y1^2*z1-q_1*y1*y4*z1-q_1*y4^2*z1-q_1*y1^2*z2-q_1*y1*y4*z2-q_1*y4^2*z2+q_1*y1*z1*z2+q_1*y4*z1*z2-q_1*y1^2*z3-q_1*y1*y4*z3-q_1*y4^2*z3+q_1*y1*z1*z3+q_1*y4*z1*z3+q_1*y1*z2*z3+q_1*y4*z2*z3-q_1*z1*z2*z3-q_1*y1^2*z4-q_1*y1*y4*z4-q_1*y4^2*z4+q_1*y1*z1*z4+q_1*y4*z1*z4+q_1*y1*z2*z4+q_1*y4*z2*z4-q_1*z1*z2*z4+q_1*y1*z3*z4+q_1*y4*z3*z4-q_1*z1*z3*z4-q_1*z2*z3*z4)*QSx([1, 4, 5, 2, 3], 'y') + (q_1*q_3*y1^2+q_1*q_3*y1*y4+q_1*q_3*y4^2+q_1*q_3*y1*y5+q_1*q_3*y4*y5+q_1*q_3*y5^2-q_1*q_3*y1*z1-q_1*q_3*y4*z1-q_1*q_3*y5*z1-q_1*q_3*y1*z2-q_1*q_3*y4*z2-q_1*q_3*y5*z2+q_1*q_3*z1*z2-q_1*q_3*y1*z3-q_1*q_3*y4*z3-q_1*q_3*y5*z3+q_1*q_3*z1*z3+q_1*q_3*z2*z3-q_1*q_3*y1*z4-q_1*q_3*y4*z4-q_1*q_3*y5*z4+q_1*q_3*z1*z4+q_1*q_3*z2*z4+q_1*q_3*z3*z4)*QSx([1, 5, 2, 3, 4], 'y') + (q_1*y1^3+q_1*y1^2*y3+q_1*y1^2*y4+q_1*y1*y3*y4+q_1*y1*y4^2+q_1*y3*y4^2+q_1*y1^2*y5+q_1*y1*y3*y5+q_1*y1*y4*y5+q_1*y3*y4*y5+q_1*y1*y5^2+q_1*y3*y5^2-2*q_1*y1^2*z1-q_1*y1*y3*z1-2*q_1*y1*y4*z1-q_1*y3*y4*z1-q_1*y4^2*z1-2*q_1*y1*y5*z1-q_1*y3*y5*z1-q_1*y4*y5*z1-q_1*y5^2*z1+q_1*y1*z1^2+q_1*y4*z1^2+q_1*y5*z1^2-2*q_1*y1^2*z2-q_1*y1*y3*z2-2*q_1*y1*y4*z2-q_1*y3*y4*z2-q_1*y4^2*z2-2*q_1*y1*y5*z2-q_1*y3*y5*z2-q_1*y4*y5*z2-q_1*y5^2*z2+3*q_1*y1*z1*z2+q_1*y3*z1*z2+2*q_1*y4*z1*z2+2*q_1*y5*z1*z2-q_1*z1^2*z2+q_1*y1*z2^2+q_1*y4*z2^2+q_1*y5*z2^2-q_1*z1*z2^2-q_1*y1^2*z3-q_1*y1*y3*z3-q_1*y1*y4*z3-q_1*y3*y4*z3-q_1*y1*y5*z3-q_1*y3*y5*z3+2*q_1*y1*z1*z3+q_1*y3*z1*z3+q_1*y4*z1*z3+q_1*y5*z1*z3-q_1*z1^2*z3+2*q_1*y1*z2*z3+q_1*y3*z2*z3+q_1*y4*z2*z3+q_1*y5*z2*z3-2*q_1*z1*z2*z3-q_1*z2^2*z3-q_1*y1^2*z4-q_1*y1*y3*z4-q_1*y1*y4*z4-q_1*y3*y4*z4-q_1*y1*y5*z4-q_1*y3*y5*z4+2*q_1*y1*z1*z4+q_1*y3*z1*z4+q_1*y4*z1*z4+q_1*y5*z1*z4-q_1*z1^2*z4+2*q_1*y1*z2*z4+q_1*y3*z2*z4+q_1*y4*z2*z4+q_1*y5*z2*z4-2*q_1*z1*z2*z4-q_1*z2^2*z4+q_1*y1*z3*z4+q_1*y3*z3*z4-q_1*z1*z3*z4-q_1*z2*z3*z4)*QSx([1, 5, 3, 2, 4], 'y') + (q_1*y1^2+q_1*y1*y4+q_1*y4^2+q_1*y1*y5+q_1*y4*y5+q_1*y5^2-q_1*y1*z1-q_1*y4*z1-q_1*y5*z1-q_1*y1*z2-q_1*y4*z2-q_1*y5*z2+q_1*z1*z2-q_1*y1*z3-q_1*y4*z3-q_1*y5*z3+q_1*z1*z3+q_1*z2*z3-q_1*y1*z4-q_1*y4*z4-q_1*y5*z4+q_1*z1*z4+q_1*z2*z4+q_1*z3*z4)*QSx([1, 5, 4, 2, 3], 'y') + (q_1*q_3*y1+q_1*q_3*y4+q_1*q_3*y5+q_1*q_3*y6-q_1*q_3*z1-q_1*q_3*z2-q_1*q_3*z3-q_1*q_3*z4)*QSx([1, 6, 2, 3, 4, 5], 'y') + (q_1*y1^2+q_1*y1*y3+q_1*y1*y4+q_1*y3*y4+q_1*y1*y5+q_1*y3*y5+q_1*y1*y6+q_1*y3*y6-2*q_1*y1*z1-q_1*y3*z1-q_1*y4*z1-q_1*y5*z1-q_1*y6*z1+q_1*z1^2-2*q_1*y1*z2-q_1*y3*z2-q_1*y4*z2-q_1*y5*z2-q_1*y6*z2+2*q_1*z1*z2+q_1*z2^2-q_1*y1*z3-q_1*y3*z3+q_1*z1*z3+q_1*z2*z3-q_1*y1*z4-q_1*y3*z4+q_1*z1*z4+q_1*z2*z4)*QSx([1, 6, 3, 2, 4, 5], 'y') + (q_1*y1+q_1*y4+q_1*y5+q_1*y6-q_1*z1-q_1*z2-q_1*z3-q_1*z4)*QSx([1, 6, 4, 2, 3, 5], 'y') + q_1*q_3*QSx([1, 7, 2, 3, 4, 5, 6], 'y') + (q_1*y1+q_1*y3-q_1*z1-q_1*z2)*QSx([1, 7, 3, 2, 4, 5, 6], 'y') + q_1*QSx([1, 7, 4, 2, 3, 5, 6], 'y') + (q_1*q_2*q_3*y1^2+q_1*q_2*q_3*y1*y2+q_1*q_2*q_3*y2^2+q_1*q_2*q_3*y1*y4+q_1*q_2*q_3*y2*y4+q_1*q_2*q_3*y4^2-q_1*q_2*q_3*y1*z1-q_1*q_2*q_3*y2*z1-q_1*q_2*q_3*y4*z1-q_1*q_2*q_3*y1*z2-q_1*q_2*q_3*y2*z2-q_1*q_2*q_3*y4*z2+q_1*q_2*q_3*z1*z2-q_1*q_2*q_3*y1*z3-q_1*q_2*q_3*y2*z3-q_1*q_2*q_3*y4*z3+q_1*q_2*q_3*z1*z3+q_1*q_2*q_3*z2*z3-q_1*q_2*q_3*y1*z4-q_1*q_2*q_3*y2*z4-q_1*q_2*q_3*y4*z4+q_1*q_2*q_3*z1*z4+q_1*q_2*q_3*z2*z4+q_1*q_2*q_3*z3*z4)*QSx([2, 1], 'y') + (q_1*q_2*q_3*y1+q_1*q_2*q_3*y2+q_1*q_2*q_3*y4+q_1*q_2*q_3*y5-q_1*q_2*q_3*z1-q_1*q_2*q_3*z2-q_1*q_2*q_3*z3-q_1*q_2*q_3*z4)*QSx([2, 1, 3, 5, 4], 'y') + q_1*q_2*q_3*QSx([2, 1, 3, 6, 4, 5], 'y') + (q_1*q_3*y1^2+q_1*q_3*y1*y2+q_1*q_3*y2^2+q_1*q_3*y1*y4+q_1*q_3*y2*y4+q_1*q_3*y4^2-q_1*q_3*y1*z1-q_1*q_3*y2*z1-q_1*q_3*y4*z1-q_1*q_3*y1*z2-q_1*q_3*y2*z2-q_1*q_3*y4*z2+q_1*q_3*z1*z2-q_1*q_3*y1*z3-q_1*q_3*y2*z3-q_1*q_3*y4*z3+q_1*q_3*z1*z3+q_1*q_3*z2*z3-q_1*q_3*y1*z4-q_1*q_3*y2*z4-q_1*q_3*y4*z4+q_1*q_3*z1*z4+q_1*q_3*z2*z4+q_1*q_3*z3*z4)*QSx([2, 4, 1, 3], 'y') + (q_1*y1^2*y3+q_1*y1*y2*y3+q_1*y2^2*y3+q_1*y1^2*y4+q_1*y1*y2*y4+q_1*y2^2*y4+q_1*y1*y3*y4+q_1*y2*y3*y4+q_1*y1*y4^2+q_1*y2*y4^2+q_1*y3*y4^2+q_1*y4^3-q_1*y1^2*z1-q_1*y1*y2*z1-q_1*y2^2*z1-q_1*y1*y3*z1-q_1*y2*y3*z1-2*q_1*y1*y4*z1-2*q_1*y2*y4*z1-q_1*y3*y4*z1-2*q_1*y4^2*z1+q_1*y1*z1^2+q_1*y2*z1^2+q_1*y4*z1^2-q_1*y1^2*z2-q_1*y1*y2*z2-q_1*y2^2*z2-q_1*y1*y3*z2-q_1*y2*y3*z2-2*q_1*y1*y4*z2-2*q_1*y2*y4*z2-q_1*y3*y4*z2-2*q_1*y4^2*z2+2*q_1*y1*z1*z2+2*q_1*y2*z1*z2+q_1*y3*z1*z2+3*q_1*y4*z1*z2-q_1*z1^2*z2+q_1*y1*z2^2+q_1*y2*z2^2+q_1*y4*z2^2-q_1*z1*z2^2-q_1*y1*y3*z3-q_1*y2*y3*z3-q_1*y1*y4*z3-q_1*y2*y4*z3-q_1*y3*y4*z3-q_1*y4^2*z3+q_1*y1*z1*z3+q_1*y2*z1*z3+q_1*y3*z1*z3+2*q_1*y4*z1*z3-q_1*z1^2*z3+q_1*y1*z2*z3+q_1*y2*z2*z3+q_1*y3*z2*z3+2*q_1*y4*z2*z3-2*q_1*z1*z2*z3-q_1*z2^2*z3-q_1*y1*y3*z4-q_1*y2*y3*z4-q_1*y1*y4*z4-q_1*y2*y4*z4-q_1*y3*y4*z4-q_1*y4^2*z4+q_1*y1*z1*z4+q_1*y2*z1*z4+q_1*y3*z1*z4+2*q_1*y4*z1*z4-q_1*z1^2*z4+q_1*y1*z2*z4+q_1*y2*z2*z4+q_1*y3*z2*z4+2*q_1*y4*z2*z4-2*q_1*z1*z2*z4-q_1*z2^2*z4+q_1*y3*z3*z4+q_1*y4*z3*z4-q_1*z1*z3*z4-q_1*z2*z3*z4)*QSx([2, 4, 3, 1], 'y') + (q_1*y1^2+q_1*y1*y2+q_1*y2^2+q_1*y1*y4+q_1*y2*y4+q_1*y4^2-q_1*y1*z1-q_1*y2*z1-q_1*y4*z1-q_1*y1*z2-q_1*y2*z2-q_1*y4*z2+q_1*z1*z2-q_1*y1*z3-q_1*y2*z3-q_1*y4*z3+q_1*z1*z3+q_1*z2*z3-q_1*y1*z4-q_1*y2*z4-q_1*y4*z4+q_1*z1*z4+q_1*z2*z4+q_1*z3*z4)*QSx([2, 4, 5, 1, 3], 'y') + (q_1*q_3*y1+q_1*q_3*y2+q_1*q_3*y4+q_1*q_3*y5-q_1*q_3*z1-q_1*q_3*z2-q_1*q_3*z3-q_1*q_3*z4)*QSx([2, 5, 1, 3, 4], 'y') + (q_1*y1^2+q_1*y1*y2+q_1*y2^2+q_1*y1*y3+q_1*y2*y3+q_1*y1*y4+q_1*y2*y4+q_1*y3*y4+q_1*y4^2+q_1*y1*y5+q_1*y2*y5+q_1*y3*y5+q_1*y4*y5+q_1*y5^2-2*q_1*y1*z1-2*q_1*y2*z1-q_1*y3*z1-2*q_1*y4*z1-2*q_1*y5*z1+q_1*z1^2-2*q_1*y1*z2-2*q_1*y2*z2-q_1*y3*z2-2*q_1*y4*z2-2*q_1*y5*z2+3*q_1*z1*z2+q_1*z2^2-q_1*y1*z3-q_1*y2*z3-q_1*y3*z3-q_1*y4*z3-q_1*y5*z3+2*q_1*z1*z3+2*q_1*z2*z3-q_1*y1*z4-q_1*y2*z4-q_1*y3*z4-q_1*y4*z4-q_1*y5*z4+2*q_1*z1*z4+2*q_1*z2*z4+q_1*z3*z4)*QSx([2, 5, 3, 1, 4], 'y') + (q_1*y1+q_1*y2+q_1*y4+q_1*y5-q_1*z1-q_1*z2-q_1*z3-q_1*z4)*QSx([2, 5, 4, 1, 3], 'y') + q_1*q_3*QSx([2, 6, 1, 3, 4, 5], 'y') + (q_1*y1+q_1*y2+q_1*y3+q_1*y4+q_1*y5+q_1*y6-2*q_1*z1-2*q_1*z2-q_1*z3-q_1*z4)*QSx([2, 6, 3, 1, 4, 5], 'y') + q_1*QSx([2, 6, 4, 1, 3, 5], 'y') + q_1*QSx([2, 7, 3, 1, 4, 5, 6], 'y') + (q_1*q_2*q_3*y1+q_1*q_2*q_3*y2+q_1*q_2*q_3*y3+q_1*q_2*q_3*y4-q_1*q_2*q_3*z1-q_1*q_2*q_3*z2-q_1*q_2*q_3*z3-q_1*q_2*q_3*z4)*QSx([3, 1, 2], 'y') + q_1*q_2*q_3*QSx([3, 1, 2, 5, 4], 'y') + (q_1*y1^2*y3+q_1*y1*y3^2+q_1*y1^2*y4+2*q_1*y1*y3*y4+q_1*y3^2*y4+q_1*y1*y4^2+q_1*y3*y4^2-q_1*y1^2*z1-2*q_1*y1*y3*z1-q_1*y3^2*z1-2*q_1*y1*y4*z1-2*q_1*y3*y4*z1-q_1*y4^2*z1+q_1*y1*z1^2+q_1*y3*z1^2+q_1*y4*z1^2-q_1*y1^2*z2-2*q_1*y1*y3*z2-q_1*y3^2*z2-2*q_1*y1*y4*z2-2*q_1*y3*y4*z2-q_1*y4^2*z2+2*q_1*y1*z1*z2+2*q_1*y3*z1*z2+2*q_1*y4*z1*z2-q_1*z1^2*z2+q_1*y1*z2^2+q_1*y3*z2^2+q_1*y4*z2^2-q_1*z1*z2^2-q_1*y1*y3*z3-q_1*y1*y4*z3-q_1*y3*y4*z3+q_1*y1*z1*z3+q_1*y3*z1*z3+q_1*y4*z1*z3-q_1*z1^2*z3+q_1*y1*z2*z3+q_1*y3*z2*z3+q_1*y4*z2*z3-q_1*z1*z2*z3-q_1*z2^2*z3-q_1*y1*y3*z4-q_1*y1*y4*z4-q_1*y3*y4*z4+q_1*y1*z1*z4+q_1*y3*z1*z4+q_1*y4*z1*z4-q_1*z1^2*z4+q_1*y1*z2*z4+q_1*y3*z2*z4+q_1*y4*z2*z4-q_1*z1*z2*z4-q_1*z2^2*z4+q_1*q_3*y1+q_1*q_3*y2+q_1*q_3*y3+q_1*q_3*y4-q_1*q_3*z1-q_1*q_3*z2-q_1*q_3*z3-q_1*q_3*z4)*QSx([3, 4, 1, 2], 'y') + (q_1*y1*y3+q_1*y2*y3+q_1*y3^2+q_1*y1*y4+q_1*y2*y4+2*q_1*y3*y4+q_1*y4^2-q_1*y1*z1-q_1*y2*z1-2*q_1*y3*z1-2*q_1*y4*z1+q_1*z1^2-q_1*y1*z2-q_1*y2*z2-2*q_1*y3*z2-2*q_1*y4*z2+2*q_1*z1*z2+q_1*z2^2-q_1*y3*z3-q_1*y4*z3+q_1*z1*z3+q_1*z2*z3-q_1*y3*z4-q_1*y4*z4+q_1*z1*z4+q_1*z2*z4)*QSx([3, 4, 2, 1], 'y') + (q_1*y1+q_1*y2+q_1*y3+q_1*y4-q_1*z1-q_1*z2-q_1*z3-q_1*z4)*QSx([3, 4, 5, 1, 2], 'y') + (q_1*y1^2+2*q_1*y1*y3+q_1*y3^2+q_1*y1*y4+q_1*y3*y4+q_1*y1*y5+q_1*y3*y5-2*q_1*y1*z1-2*q_1*y3*z1-q_1*y4*z1-q_1*y5*z1+q_1*z1^2-2*q_1*y1*z2-2*q_1*y3*z2-q_1*y4*z2-q_1*y5*z2+2*q_1*z1*z2+q_1*z2^2-q_1*y1*z3-q_1*y3*z3+q_1*z1*z3+q_1*z2*z3-q_1*y1*z4-q_1*y3*z4+q_1*z1*z4+q_1*z2*z4+q_1*q_3)*QSx([3, 5, 1, 2, 4], 'y') + (q_1*y1+q_1*y2+2*q_1*y3+q_1*y4+q_1*y5-2*q_1*z1-2*q_1*z2-q_1*z3-q_1*z4)*QSx([3, 5, 2, 1, 4], 'y') + q_1*QSx([3, 5, 4, 1, 2], 'y') + (q_1*y1+q_1*y3-q_1*z1-q_1*z2)*QSx([3, 6, 1, 2, 4, 5], 'y') + q_1*QSx([3, 6, 2, 1, 4, 5], 'y') + (q_3*y4^4-q_3*y4^3*z1-q_3*y4^3*z2+q_3*y4^2*z1*z2-q_3*y4^3*z3+q_3*y4^2*z1*z3+q_3*y4^2*z2*z3-q_3*y4*z1*z2*z3-q_3*y4^3*z4+q_3*y4^2*z1*z4+q_3*y4^2*z2*z4-q_3*y4*z1*z2*z4+q_3*y4^2*z3*z4-q_3*y4*z1*z3*z4-q_3*y4*z2*z3*z4+q_3*z1*z2*z3*z4)*QSx([4, 1, 2, 3], 'y') + (y1*y4^4+y3*y4^4-y1*y4^3*z1-y3*y4^3*z1-y4^4*z1+y4^3*z1^2-y1*y4^3*z2-y3*y4^3*z2-y4^4*z2+y1*y4^2*z1*z2+y3*y4^2*z1*z2+2*y4^3*z1*z2-y4^2*z1^2*z2+y4^3*z2^2-y4^2*z1*z2^2-y1*y4^3*z3-y3*y4^3*z3+y1*y4^2*z1*z3+y3*y4^2*z1*z3+y4^3*z1*z3-y4^2*z1^2*z3+y1*y4^2*z2*z3+y3*y4^2*z2*z3+y4^3*z2*z3-y1*y4*z1*z2*z3-y3*y4*z1*z2*z3-2*y4^2*z1*z2*z3+y4*z1^2*z2*z3-y4^2*z2^2*z3+y4*z1*z2^2*z3-y1*y4^3*z4-y3*y4^3*z4+y1*y4^2*z1*z4+y3*y4^2*z1*z4+y4^3*z1*z4-y4^2*z1^2*z4+y1*y4^2*z2*z4+y3*y4^2*z2*z4+y4^3*z2*z4-y1*y4*z1*z2*z4-y3*y4*z1*z2*z4-2*y4^2*z1*z2*z4+y4*z1^2*z2*z4-y4^2*z2^2*z4+y4*z1*z2^2*z4+y1*y4^2*z3*z4+y3*y4^2*z3*z4-y1*y4*z1*z3*z4-y3*y4*z1*z3*z4-y4^2*z1*z3*z4+y4*z1^2*z3*z4-y1*y4*z2*z3*z4-y3*y4*z2*z3*z4-y4^2*z2*z3*z4+y1*z1*z2*z3*z4+y3*z1*z2*z3*z4+2*y4*z1*z2*z3*z4-z1^2*z2*z3*z4+y4*z2^2*z3*z4-z1*z2^2*z3*z4)*QSx([4, 1, 3, 2], 'y') + (y4^4-y4^3*z1-y4^3*z2+y4^2*z1*z2-y4^3*z3+y4^2*z1*z3+y4^2*z2*z3-y4*z1*z2*z3-y4^3*z4+y4^2*z1*z4+y4^2*z2*z4-y4*z1*z2*z4+y4^2*z3*z4-y4*z1*z3*z4-y4*z2*z3*z4+z1*z2*z3*z4)*QSx([4, 1, 5, 2, 3], 'y') + (y4^4-y4^3*z1-y4^3*z2+y4^2*z1*z2-y4^3*z3+y4^2*z1*z3+y4^2*z2*z3-y4*z1*z2*z3-y4^3*z4+y4^2*z1*z4+y4^2*z2*z4-y4*z1*z2*z4+y4^2*z3*z4-y4*z1*z3*z4-y4*z2*z3*z4+z1*z2*z3*z4)*QSx([4, 2, 3, 1], 'y') + (q_1*y1+q_1*y3+q_1*y4+q_1*y5-q_1*z1-q_1*z2-q_1*z3-q_1*z4)*QSx([4, 5, 1, 2, 3], 'y') + q_1*QSx([4, 5, 2, 1, 3], 'y') + q_1*QSx([4, 6, 1, 2, 3, 5], 'y') + (q_3*y4^3+q_3*y4^2*y5+q_3*y4*y5^2+q_3*y5^3-q_3*y4^2*z1-q_3*y4*y5*z1-q_3*y5^2*z1-q_3*y4^2*z2-q_3*y4*y5*z2-q_3*y5^2*z2+q_3*y4*z1*z2+q_3*y5*z1*z2-q_3*y4^2*z3-q_3*y4*y5*z3-q_3*y5^2*z3+q_3*y4*z1*z3+q_3*y5*z1*z3+q_3*y4*z2*z3+q_3*y5*z2*z3-q_3*z1*z2*z3-q_3*y4^2*z4-q_3*y4*y5*z4-q_3*y5^2*z4+q_3*y4*z1*z4+q_3*y5*z1*z4+q_3*y4*z2*z4+q_3*y5*z2*z4-q_3*z1*z2*z4+q_3*y4*z3*z4+q_3*y5*z3*z4-q_3*z1*z3*z4-q_3*z2*z3*z4)*QSx([5, 1, 2, 3, 4], 'y') + (y1*y4^3+y3*y4^3+y1*y4^2*y5+y3*y4^2*y5+y1*y4*y5^2+y3*y4*y5^2+y1*y5^3+y3*y5^3-y1*y4^2*z1-y3*y4^2*z1-y4^3*z1-y1*y4*y5*z1-y3*y4*y5*z1-y4^2*y5*z1-y1*y5^2*z1-y3*y5^2*z1-y4*y5^2*z1-y5^3*z1+y4^2*z1^2+y4*y5*z1^2+y5^2*z1^2-y1*y4^2*z2-y3*y4^2*z2-y4^3*z2-y1*y4*y5*z2-y3*y4*y5*z2-y4^2*y5*z2-y1*y5^2*z2-y3*y5^2*z2-y4*y5^2*z2-y5^3*z2+y1*y4*z1*z2+y3*y4*z1*z2+2*y4^2*z1*z2+y1*y5*z1*z2+y3*y5*z1*z2+2*y4*y5*z1*z2+2*y5^2*z1*z2-y4*z1^2*z2-y5*z1^2*z2+y4^2*z2^2+y4*y5*z2^2+y5^2*z2^2-y4*z1*z2^2-y5*z1*z2^2-y1*y4^2*z3-y3*y4^2*z3-y1*y4*y5*z3-y3*y4*y5*z3-y1*y5^2*z3-y3*y5^2*z3+y1*y4*z1*z3+y3*y4*z1*z3+y4^2*z1*z3+y1*y5*z1*z3+y3*y5*z1*z3+y4*y5*z1*z3+y5^2*z1*z3-y4*z1^2*z3-y5*z1^2*z3+y1*y4*z2*z3+y3*y4*z2*z3+y4^2*z2*z3+y1*y5*z2*z3+y3*y5*z2*z3+y4*y5*z2*z3+y5^2*z2*z3-y1*z1*z2*z3-y3*z1*z2*z3-2*y4*z1*z2*z3-2*y5*z1*z2*z3+z1^2*z2*z3-y4*z2^2*z3-y5*z2^2*z3+z1*z2^2*z3-y1*y4^2*z4-y3*y4^2*z4-y1*y4*y5*z4-y3*y4*y5*z4-y1*y5^2*z4-y3*y5^2*z4+y1*y4*z1*z4+y3*y4*z1*z4+y4^2*z1*z4+y1*y5*z1*z4+y3*y5*z1*z4+y4*y5*z1*z4+y5^2*z1*z4-y4*z1^2*z4-y5*z1^2*z4+y1*y4*z2*z4+y3*y4*z2*z4+y4^2*z2*z4+y1*y5*z2*z4+y3*y5*z2*z4+y4*y5*z2*z4+y5^2*z2*z4-y1*z1*z2*z4-y3*z1*z2*z4-2*y4*z1*z2*z4-2*y5*z1*z2*z4+z1^2*z2*z4-y4*z2^2*z4-y5*z2^2*z4+z1*z2^2*z4+y1*y4*z3*z4+y3*y4*z3*z4+y1*y5*z3*z4+y3*y5*z3*z4-y1*z1*z3*z4-y3*z1*z3*z4-y4*z1*z3*z4-y5*z1*z3*z4+z1^2*z3*z4-y1*z2*z3*z4-y3*z2*z3*z4-y4*z2*z3*z4-y5*z2*z3*z4+2*z1*z2*z3*z4+z2^2*z3*z4)*QSx([5, 1, 3, 2, 4], 'y') + (y4^3+y4^2*y5+y4*y5^2+y5^3-y4^2*z1-y4*y5*z1-y5^2*z1-y4^2*z2-y4*y5*z2-y5^2*z2+y4*z1*z2+y5*z1*z2-y4^2*z3-y4*y5*z3-y5^2*z3+y4*z1*z3+y5*z1*z3+y4*z2*z3+y5*z2*z3-z1*z2*z3-y4^2*z4-y4*y5*z4-y5^2*z4+y4*z1*z4+y5*z1*z4+y4*z2*z4+y5*z2*z4-z1*z2*z4+y4*z3*z4+y5*z3*z4-z1*z3*z4-z2*z3*z4)*QSx([5, 1, 4, 2, 3], 'y') + (y4^3+y4^2*y5+y4*y5^2+y5^3-y4^2*z1-y4*y5*z1-y5^2*z1-y4^2*z2-y4*y5*z2-y5^2*z2+y4*z1*z2+y5*z1*z2-y4^2*z3-y4*y5*z3-y5^2*z3+y4*z1*z3+y5*z1*z3+y4*z2*z3+y5*z2*z3-z1*z2*z3-y4^2*z4-y4*y5*z4-y5^2*z4+y4*z1*z4+y5*z1*z4+y4*z2*z4+y5*z2*z4-z1*z2*z4+y4*z3*z4+y5*z3*z4-z1*z3*z4-z2*z3*z4)*QSx([5, 2, 3, 1, 4], 'y') + (q_3*y4^2+q_3*y4*y5+q_3*y5^2+q_3*y4*y6+q_3*y5*y6+q_3*y6^2-q_3*y4*z1-q_3*y5*z1-q_3*y6*z1-q_3*y4*z2-q_3*y5*z2-q_3*y6*z2+q_3*z1*z2-q_3*y4*z3-q_3*y5*z3-q_3*y6*z3+q_3*z1*z3+q_3*z2*z3-q_3*y4*z4-q_3*y5*z4-q_3*y6*z4+q_3*z1*z4+q_3*z2*z4+q_3*z3*z4)*QSx([6, 1, 2, 3, 4, 5], 'y') + (y1*y4^2+y3*y4^2+y1*y4*y5+y3*y4*y5+y1*y5^2+y3*y5^2+y1*y4*y6+y3*y4*y6+y1*y5*y6+y3*y5*y6+y1*y6^2+y3*y6^2-y1*y4*z1-y3*y4*z1-y4^2*z1-y1*y5*z1-y3*y5*z1-y4*y5*z1-y5^2*z1-y1*y6*z1-y3*y6*z1-y4*y6*z1-y5*y6*z1-y6^2*z1+y4*z1^2+y5*z1^2+y6*z1^2-y1*y4*z2-y3*y4*z2-y4^2*z2-y1*y5*z2-y3*y5*z2-y4*y5*z2-y5^2*z2-y1*y6*z2-y3*y6*z2-y4*y6*z2-y5*y6*z2-y6^2*z2+y1*z1*z2+y3*z1*z2+2*y4*z1*z2+2*y5*z1*z2+2*y6*z1*z2-z1^2*z2+y4*z2^2+y5*z2^2+y6*z2^2-z1*z2^2-y1*y4*z3-y3*y4*z3-y1*y5*z3-y3*y5*z3-y1*y6*z3-y3*y6*z3+y1*z1*z3+y3*z1*z3+y4*z1*z3+y5*z1*z3+y6*z1*z3-z1^2*z3+y1*z2*z3+y3*z2*z3+y4*z2*z3+y5*z2*z3+y6*z2*z3-2*z1*z2*z3-z2^2*z3-y1*y4*z4-y3*y4*z4-y1*y5*z4-y3*y5*z4-y1*y6*z4-y3*y6*z4+y1*z1*z4+y3*z1*z4+y4*z1*z4+y5*z1*z4+y6*z1*z4-z1^2*z4+y1*z2*z4+y3*z2*z4+y4*z2*z4+y5*z2*z4+y6*z2*z4-2*z1*z2*z4-z2^2*z4+y1*z3*z4+y3*z3*z4-z1*z3*z4-z2*z3*z4)*QSx([6, 1, 3, 2, 4, 5], 'y') + (y4^2+y4*y5+y5^2+y4*y6+y5*y6+y6^2-y4*z1-y5*z1-y6*z1-y4*z2-y5*z2-y6*z2+z1*z2-y4*z3-y5*z3-y6*z3+z1*z3+z2*z3-y4*z4-y5*z4-y6*z4+z1*z4+z2*z4+z3*z4)*QSx([6, 1, 4, 2, 3, 5], 'y') + (y4^2+y4*y5+y5^2+y4*y6+y5*y6+y6^2-y4*z1-y5*z1-y6*z1-y4*z2-y5*z2-y6*z2+z1*z2-y4*z3-y5*z3-y6*z3+z1*z3+z2*z3-y4*z4-y5*z4-y6*z4+z1*z4+z2*z4+z3*z4)*QSx([6, 2, 3, 1, 4, 5], 'y') + (q_3*y4+q_3*y5+q_3*y6+q_3*y7-q_3*z1-q_3*z2-q_3*z3-q_3*z4)*QSx([7, 1, 2, 3, 4, 5, 6], 'y') + (y1*y4+y3*y4+y1*y5+y3*y5+y1*y6+y3*y6+y1*y7+y3*y7-y1*z1-y3*z1-y4*z1-y5*z1-y6*z1-y7*z1+z1^2-y1*z2-y3*z2-y4*z2-y5*z2-y6*z2-y7*z2+2*z1*z2+z2^2-y1*z3-y3*z3+z1*z3+z2*z3-y1*z4-y3*z4+z1*z4+z2*z4)*QSx([7, 1, 3, 2, 4, 5, 6], 'y') + (y4+y5+y6+y7-z1-z2-z3-z4)*QSx([7, 1, 4, 2, 3, 5, 6], 'y') + (y4+y5+y6+y7-z1-z2-z3-z4)*QSx([7, 2, 3, 1, 4, 5, 6], 'y') + q_3*QSx([8, 1, 2, 3, 4, 5, 6, 7], 'y') + (y1+y3-z1-z2)*QSx([8, 1, 3, 2, 4, 5, 6, 7], 'y') + QSx([8, 1, 4, 2, 3, 5, 6, 7], 'y') + QSx([8, 2, 3, 1, 4, 5, 6, 7], 'y')
```
This output of roughly 17,000 characters took about 2 seconds to compute on my laptop. **Again note that quantum double computations are technically conjectural, but a proof is likely forthcoming.**

[Homepage of schubmult](http://schubmult.org/)