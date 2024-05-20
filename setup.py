from setuptools import setup, find_packages

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="schubmult",
    version="1.3.8",
    description="Computing Littlewood-Richardson coefficients of Schubert polynomials",
	long_description=long_description,
	long_description_content_type='text/markdown',
    url="https://github.com/matthematics/schubmult",
    author="Matt Samuel",
    author_email="schubmult@gmail.com",
    license="GNU",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
		"PuLP>=2.7.0",
		"symengine>=0.11.0",
		"sympy>=1.12",
		"psutil",
		"cachetools",
		"sortedcontainers"
    ],
    entry_points={"console_scripts": ["schubmult_py=schubmult.schubmult_py.__main__:main",
	"schubmult_double=schubmult.schubmult_double.__main__:main",
	"schubmult_yz=schubmult.schubmult_yz.__main__:main",
	"schubmult_q=schubmult.schubmult_q.__main__:main",
	"schubmult_q_double=schubmult.schubmult_q_double.__main__:main",
	"schubmult_q_yz=schubmult.schubmult_q_yz.__main__:main"
	]},
)