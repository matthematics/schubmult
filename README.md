# schubmult
Program for computing Littlewood-Richardson coefficients of Schubert polynomials

This is a python script written by Matt Samuel for computing Littlewood-Richardson coefficients of Schubert polynomials. It has the same command line syntax as the program "schubmult" in lrcalc by Anders Buch. The speed of the code is comparable and sometimes better than Buch's schubmult as of the time of this writing, September 10th, 2023. Example:

python3 schubmult_py.py 1 2 4 9 11 6 8 12 3 5 7 10 - 6 8 1 2 3 4 7 10 12 14 5 9 11 13

or use the executable if it is installed. This has a runtime of roughly 1 minute on my machine, whereas schubmult takes 1 minute, 19 seconds. Runtime will vary tremendously by case. This problem is #P-hard. Though the result is always nonnegative and the problem is in GapP, it is not known to be in #P at this time.

schubmult_yz.py is for multiplying double Schubert polynomials in different sets of coefficient variables (labeled y and z), and schubmult_double.py is for multiplying double Schubert polynomials in the same set of coefficient variables. Both have the same command line syntax as schubmult.py. schubmult_double.py displays the result with nonnegative coefficients in terms of the negative simple roots. Both are of course slower than schubmult.py, and expressing the result positively for schubmult_double.py slows it down even more.

Both schubmult_yz.py and schubmult_double.py require symengine. For versions that are much slower but use only sympy, use the _sympy scripts.
