from distutils.core import setup
import schubmult.schubmult_py
import schubmult.schubmult_yz
import schubmult.schubmult_double


setup(
    name="schubmult",
    version="1.0.0",
    description="Computing Littlewood-Richardson coefficients of Schubert polynomials",
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