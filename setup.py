from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="schubmult",
    version="1.0.12",
    description="Computing Littlewood-Richardson coefficients of Schubert polynomials",
	long_description=long_description,
	long_description_content_type='text/markdown',
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