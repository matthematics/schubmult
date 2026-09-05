import subprocess
import sys

import pytest

# Every ring submodule must be importable as the very first schubmult import.
# Regression guard: schubert_ring <-> quantum_schubert_ring <-> parabolic_quantum_double_schubert_ring
# formed an import cycle, so entering it at the parabolic module raised ImportError.
RING_MODULES = [
    "schubmult.rings.schubert.base_schubert_ring",
    "schubmult.rings.schubert.double_schubert_ring",
    "schubmult.rings.schubert.parabolic_quantum_double_schubert_ring",
    "schubmult.rings.schubert.parabolic_quantum_schubert_ring",
    "schubmult.rings.schubert.quantum_double_schubert_ring",
    "schubmult.rings.schubert.quantum_schubert_ring",
    "schubmult.rings.schubert.schubert_ring",
]


@pytest.mark.parametrize("module", RING_MODULES)
def test_module_imports_in_isolation(module):
    result = subprocess.run(
        [sys.executable, "-c", f"import {module}"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr


PARABOLIC_NAMES = [
    "ParabolicQuantumDoubleSchubertElement",
    "ParabolicQuantumDoubleSchubertRing",
    "QPDSx",
    "QPSx",
    "make_parabolic_quantum_basis",
]


@pytest.mark.parametrize("name", PARABOLIC_NAMES)
def test_parabolic_name_exported(name):
    import schubmult.rings.schubert as rings

    assert getattr(rings, name) is not None


def test_quantum_cohomology_projective_space():
    """QH^*(P^3): h^4 = q_1 and h^5 = q_1 h."""
    from schubmult.combinatorics.permutation import Permutation
    from schubmult.rings.schubert import QPSx
    from schubmult.symbolic.poly.variables import GeneratingSet

    q_var = GeneratingSet("q")
    ring = QPSx(1, 3)
    h = ring([2, 1, 3, 4])

    assert dict((h**4).kill_ideal().items()) == {Permutation([]): q_var[1]}
    assert dict((h**5).kill_ideal().items()) == {Permutation([2, 1]): q_var[1]}


def test_quantum_cohomology_grassmannian_2_4():
    """QH^*(Gr(2,4)): sigma_1^2 = sigma_2 + sigma_11, sigma_1^4 = 2 sigma_22 + 2 q_1."""
    from schubmult.combinatorics.permutation import Permutation
    from schubmult.rings.schubert import QPSx
    from schubmult.symbolic.poly.variables import GeneratingSet

    q_var = GeneratingSet("q")
    ring = QPSx(2, 2)
    sigma1 = ring([1, 3, 2, 4])

    assert dict((sigma1**2).kill_ideal().items()) == {Permutation([1, 4, 2, 3]): 1, Permutation([2, 3, 1]): 1}
    assert dict((sigma1**4).kill_ideal().items()) == {Permutation([3, 4, 1, 2]): 2, Permutation([]): 2 * q_var[1]}


def test_single_and_double_parabolic_constructors():
    from schubmult.rings.schubert import make_parabolic_quantum_basis
    from schubmult.rings.schubert.parabolic_quantum_schubert_ring import make_single_parabolic_quantum_basis
    from schubmult.symbolic.poly.variables import GeneratingSet

    single = make_single_parabolic_quantum_basis((1, 3))
    double = make_parabolic_quantum_basis((1, 3), GeneratingSet("y"))

    assert single([2, 1, 3, 4]) is not None
    assert double([2, 1, 3, 4]) is not None
