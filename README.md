# schubmult

## Program and package for computing Littlewood-Richardson coefficients of Schubert polynomials

This is a set of python scripts written by Matt Samuel for computing Littlewood-Richardson coefficients of (ordinary or double) Schubert polynomials. Since version 1.3.3, it also handles (double) quantum Schubert polynomials, if double then either in the same set or different sets of coefficient variables; that is to say it compute the (equivariant/mixed) Gromov-Witten invariants of the complete flag variety. It has the same command line syntax as the program "schubmult" in lrcalc by Anders Buch. Example:

```
schubmult_py 1 2 4 9 11 6 8 12 3 5 7 10 - 6 8 1 2 3 4 7 10 12 14 5 9 11 13  
schubmult_double 1 3 4 6 2 5 - 2 1 5 7 3 4 6  
schubmult_yz 1 3 4 6 2 5 - 2 1 5 7 3 4 6 --display-positive
schubmult_q 5 1 4 3 2 - 5 1 3 4 2
schubmult_q_double 5 1 4 3 2 - 5 1 3 4 2
schubmult_q_yz 5 1 4 3 2 - 2 5 1 3 4 --display-positive
```

The same execution with the Lehmer code:

```
schubmult_py -code 0 0 1 5 6 2 3 4 - 5 6 0 0 0 0 1 2 3 4
schubmult_double -code 0 1 1 2 - 1 0 2 3
schubmult_yz -code 0 1 1 2 - 1 0 2 3 --display-positive
schubmult_q -code 4 0 2 1 - 4 0 1 1
schubmult_q_double -code 4 0 2 1 - 4 0 1 1
schubmult_q_yz -code 4 0 2 1 - 1 3 --display-positive
```

For coproducts:
```
schubmult_py -coprod 1 3 5 7 2 4 6 - 2 4
schubmult_double -coprod 1 3 5 7 2 4 6 - 2 4
schubmult_yz -coprod 1 3 5 7 2 4 6 - 2 4 --display-positive
```
or
```
schubmult_py -code -coprod 0 1 2 3 - 2 4
schubmult_double -code -coprod 0 1 2 3 - 2 4
schubmult_yz -code -coprod 0 1 2 3 - 2 4 --display-positive
```

Since version 1.3.7, schubmult_q_yz has a feature for displaying the coefficients of the divided difference operators in the evaluation of the quantum double Schubert polynomials on the commuting difference operators of Fomin, Gelfand, and Postnikov. It is necessary to cap the value of n in the group S_n we are working in because as n increases the expression does not stabilize.
```
schubmult_q_yz -nil-hecke 6 -code 2 2 --display-positive
```

Runtime will vary tremendously by case. The general problem is #P-hard. Though the result is always nonnegative (which at least is known for schubmult_py, schubmult_q, schubmult_double, and schubmult_q_double) and the problem is in GapP, it is not known to be in #P at this time.

schubmult_py is for multiplying ordinary Schubert polynomials. schubmult_yz is for multiplying double Schubert polynomials in different sets of coefficient variables (labeled y and z), and schubmult_double is for multiplying double Schubert polynomials in the same set of coefficient variables. Similarly, schubmult_q is for multiplying quantum Schubert polynomials, schubmult_q_double is for multiplying quantum double Schubert polynomials in the same set of coefficient variables, and schubmult_q_yz is for multiplying quantum double Schubert polynomials in different sets of coefficient variables, or in other words it computes the Gromov-Witten invariants, equivariant Gromov-Witten invariants, and (mixed?) equivariant Gromov-Witten invariants of the complete flag variety. All have the same command line syntax as schubmult, except when using the -code option. schubmult_double/schubmult_q_double display the result with nonnegative coefficients in terms of the negative simple roots (and the q variables), and schubmult_yz and schubmult_q_yz optionally display the result positively in terms of y_i-z_j (and q) with the --display-positive option.

New in version 1.1.0, schubmult_xx -coprod allows you to split (double) Schubert polynomials along certain indices (not available for schubmult_q). It takes one permutation as an argument, followed by a dash -, then the set of indices you would like to split on. These coefficients are always nonnegative since they occur as product coefficients (this is actually how they are computed).

When imported as a python package, the relevant packages are schubmult.perm_lib, which has various permutation manipulation functions, and three modules that have functions of the same name (function name is "schubmult"): schubmult.schubmult_py, schubmult.schubmult_yz, schubmult.schubmult_double. Function takes a permutation dictionary (keys are tuples of ints, which must be trimmed permutations, and values are either integers or symengine values, which can also be integers) as well as a permutation as its second argument, which is the (double) Schubert polynomial to multiply by. Returns a dictionary of the same form with the coefficients.

```
from schubmult.schubmult_yz import schubmult  
  
coeff_dict = schubmult({(1,3,4,6,2,5): 1},(2,1,5,7,3,4,6)) # outputs dictionary with results  
```


```
from schubmult.schubmult_py import schubmult  
  
coeff_dict = schubmult({(1,3,4,6,2,5): 1},(2,1,5,7,3,4,6))
```

Version 1.0.18 adds the command line argument --display-positive to schubmult_yz (and version 1.3.3 adds --display-positive to schubmult_q_yz), which displays the result positively (if possible, this is still only always possible conjecturally). It will fail and print out the offending case if it finds a counterexample. This is highly processor intensive.

![](https://raw.githubusercontent.com/matthematics/schubmult/main/positive_image.png)

[Homepage of schubmult](http://schubmult.org/)