# schubmult
Program for computing Littlewood-Richardson coefficients of Schubert polynomials

This is a python script written by Matt Samuel for computing Littlewood-Richardson coefficients of Schubert polynomials. It has the same command line syntax as the program "schubmult" in lrcalc by Anders Buch. The speed of the code is comparable and sometimes better than Buch's schubmult as of the time of this writing, September 10th, 2023. Example:

python3 schubmult.py 1 2 4 9 11 6 8 12 3 5 7 10 - 6 8 1 2 3 4 7 10 12 14 5 9 11 13

This has a runtime of roughly 1 minute on my machine, whereas schubmult takes 1 minute, 19 seconds. Speed will vary tremendously by case. This problem is #P-hard. Though the result is always nonnegative and the problem is in GapP, it is not known to be in #P at this time.
