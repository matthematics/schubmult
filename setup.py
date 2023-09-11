from distutils.core import setup
import schubmult.schubmult_py
import schubmult.schubmult_yz
import schubmult.schubmult_double


setup(
    name="schubmult",
    version="1.0.0",
    description="Computing Littlewood-Richardson coefficients of Schubert polynomials",
	long_description="""This package is primarily intended for installing the executables schubmult_py, schubmult_yz, and schubmult_double.
	It exposes some permutation functions as well as three versions of the 'schubmult' function, which are not currently intended for public use.""",
    url="https://github.com/matthematics/schubmult",
    author="Matt Samuel",
#    author_email="info@realpython.com",
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