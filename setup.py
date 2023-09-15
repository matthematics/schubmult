from distutils.core import setup
import schubmult.schubmult_py
import schubmult.schubmult_yz
import schubmult.schubmult_double


setup(
    name="schubmult",
    version="1.0.7",
    description="Computing Littlewood-Richardson coefficients of Schubert polynomials",
	long_description="""Program for computing Littlewood-Richardson coefficients of Schubert polynomials

This is a python package written by Matt Samuel for computing Littlewood-Richardson coefficients of Schubert polynomials. It has three scripts installed as executables with the same command line syntax as the program "schubmult" in lrcalc by Anders Buch. Example:

schubmult_py 1 2 4 9 11 6 8 12 3 5 7 10 - 6 8 1 2 3 4 7 10 12 14 5 9 11 13

The general problem is #P-hard. Though the result is always nonnegative and the problem is in GapP, it is not known to be in #P at this time.

Do not try this example with schubmult_yz or schubmult_double, which are the other two scripts, because it is too large for most machines.

schubmult_yz is for multiplying double Schubert polynomials in different sets of coefficient variables (labeled y and z), and schubmult_double is for multiplying double Schubert polynomials in the same set of coefficient variables. Both have the same command line syntax as schubmult. schubmult_double displays the result with nonnegative coefficients in terms of the negative simple roots. Both are of course slower than schubmult_py, and expressing the result positively for schubmult_double slows it down even more.

When imported as a python package, the relevant packages are schubmult.perm_lib, which has various permutation manipulation functions, and three modules that have functions of the same name (function name is "schubmult"): schubmult.schubmult_py, schubmult.schubmult_yz, schubmult.schubmult_double. Function takes a permutation dictionary (keys are tuples of ints, which must be trimmed permutations, and values are either integers or symengine values, which can also be integers) as well as a permutation as its second argument, which is the (double) Schubert polynomial to multiply by. Returns a dictionary of the same form with the coefficients.
	""",
    url="https://github.com/matthematics/schubmult",
    author="Matt Samuel",
    author_email="matthematics@gmail.com",
    license="GNU",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    packages=["schubmult","schubmult.schubmult_py","schubmult.schubmult_double","schubmult.schubmult_yz"],
    include_package_data=True,
    install_requires=[
        "symengine", "numpy"
    ],
    entry_points={"console_scripts": ["schubmult_py=schubmult.schubmult_py.__main__:main",
	"schubmult_double=schubmult.schubmult_double.__main__:main",
	"schubmult_yz=schubmult.schubmult_yz.__main__:main"
	]},
)