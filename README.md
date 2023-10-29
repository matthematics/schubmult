# schubmult

## Program and package for computing Littlewood-Richardson coefficients of Schubert polynomials

This is a set of python scripts written by Matt Samuel for computing Littlewood-Richardson coefficients of (ordinary or double) Schubert polynomials. It has the same command line syntax as the program "schubmult" in lrcalc by Anders Buch. Example:

```
schubmult_py 1 2 4 9 11 6 8 12 3 5 7 10 - 6 8 1 2 3 4 7 10 12 14 5 9 11 13  
schubmult_double 1 3 4 6 2 5 - 2 1 5 7 3 4 6  
schubmult_yz 1 3 4 6 2 5 - 2 1 5 7 3 4 6

schubmult_py -coprod 1 3 5 7 2 4 6 - 2 4
```

Runtime will vary tremendously by case. The general problem is #P-hard. Though the result is always nonnegative and the problem is in GapP, it is not known to be in #P at this time.

schubmult_py is for multiplying ordinary Schubert polynomials. schubmult_yz is for multiplying double Schubert polynomials in different sets of coefficient variables (labeled y and z), and schubmult_double is for multiplying double Schubert polynomials in the same set of coefficient variables. Both have the same command line syntax as schubmult. schubmult_double displays the result with nonnegative coefficients in terms of the negative simple roots. Both are of course slower than schubmult_py, and expressing the result positively for schubmult_double slows it down even more.

New in version 1.1.0, schubmult_py -coprod allows you to split Schubert polynomials along certain indices. It takes one permutation as an argument, followed by a dash -, then the set of indices you would like to split on. These coefficients are always nonnegative since they occur as product coefficients (this is actually how they are computed).

When imported as a python package, the relevant packages are schubmult.perm_lib, which has various permutation manipulation functions, and three modules that have functions of the same name (function name is "schubmult"): schubmult.schubmult_py, schubmult.schubmult_yz, schubmult.schubmult_double. Function takes a permutation dictionary (keys are tuples of ints, which must be trimmed permutations, and values are either integers or symengine values, which can also be integers) as well as a permutation as its second argument, which is the (double) Schubert polynomial to multiply by. Returns a dictionary of the same form with the coefficients.

```
from schubmult.schubmult_yz import schubmult  
  
coeff_dict = schubmult({(1,3,4,6,2,5): 1},(2,1,5,7,3,4,6)) # outputs dictionary with results  
```


```
from schubmult.schubmult_py import schubmult  
  
coeff_dict = schubmult({(1,3,4,6,2,5): 1},(2,1,5,7,3,4,6))
```

Note versions 1.0.15 and prior had a bug that failed to upgrade the executable. The coefficients it computed were correct, but it was not using the updated version.

Version 1.0.18 adds the command line argument --display-positive to schubmult_yz, which displays the result positively (if possible, this is still only always possible conjecturally). This is highly processor intensive.

![](https://raw.githubusercontent.com/matthematics/schubmult/main/positive_image.png)

[Homepage of schubmult](http://schubmult.org/)