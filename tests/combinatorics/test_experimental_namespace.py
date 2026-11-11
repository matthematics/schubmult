import importlib
import subprocess
import sys

import pytest

import schubmult
from schubmult.combinatorics import experimental


_MODULE_EXPORTS = {
    "chute_move_element": ("ChuteMoveElement",),
    "hecke_plactic": ("HeckePlactic",),
    "hpd": ("HPD", "HPDTile"),
    "inversions_tableau": ("InversionsTableau",),
    "quasi_crystal_graph": ("QuasiCrystalGraph", "QuasiCrystalGraphTensor"),
    "root_tableau": ("RootTableau",),
}


@pytest.mark.parametrize(("module", "names"), _MODULE_EXPORTS.items())
def test_experimental_exports(module, names):
    module_name = f"schubmult.combinatorics.experimental.{module}"
    implementation = importlib.import_module(module_name)
    for name in names:
        exported = getattr(experimental, name)
        assert exported is getattr(implementation, name)
        assert exported.__module__ == module_name


def test_experimental_exports_load_lazily():
    code = """
import sys
import schubmult.combinatorics.experimental as experimental

module_name = "schubmult.combinatorics.experimental.hpd"
assert module_name not in sys.modules
assert experimental.HPD.__module__ == module_name
assert module_name in sys.modules
"""
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("module", _MODULE_EXPORTS)
def test_old_module_paths_are_not_importable(module):
    module_name = f"schubmult.combinatorics.{module}"
    with pytest.raises(ModuleNotFoundError) as exc_info:
        importlib.import_module(module_name)
    assert exc_info.value.name == module_name


@pytest.mark.parametrize("name", [name for names in _MODULE_EXPORTS.values() for name in names])
def test_experimental_names_are_not_top_level(name):
    assert name not in schubmult.__all__
    assert name not in dir(schubmult)
    with pytest.raises(AttributeError):
        getattr(schubmult, name)


def test_chute_move_ring_uses_experimental_element():
    from schubmult import ChuteMoveRing
    from schubmult.combinatorics.experimental import ChuteMoveElement

    ring = ChuteMoveRing()
    assert isinstance(ring.zero_monom, ChuteMoveElement)
    assert ring.one * ring.one == ring.one


def test_rc_graph_eg_error_points_to_experimental_tableau():
    from schubmult import RCGraph

    with pytest.raises(NotImplementedError, match=r"schubmult\.combinatorics\.experimental\.root_tableau\.RootTableau\.from_rc_graph"):
        RCGraph([(1,)]).edelman_greene()


def test_alt_rc_graph_ring_uses_experimental_tableau():
    from schubmult import RCGraph
    from schubmult.rings.combinatorial.alt_rc_graph_ring import AltRCGraphRing

    ring = AltRCGraphRing()
    rc = RCGraph([(1,)])
    assert dict(ring.rc_single_product(rc, rc)) == {RCGraph([(2, 1)]): 1}
