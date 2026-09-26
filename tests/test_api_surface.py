"""Guards against API-surface breakage that correctness tests do not catch.

The bug fixed in 4.2.0 was an ImportError, not a wrong answer: every module
imported fine in the order the test suite happened to use, but entering the
import graph at ``parabolic_quantum_double_schubert_ring`` hit a cycle. Nothing
in the suite imported it, so nothing failed.

Two properties are checked here, both of which that bug violated:

1. Every name the package advertises can actually be retrieved. The top-level
   ``__getattr__`` imports lazily and converts any failure into AttributeError,
   so a broken module stays invisible until someone touches the name.
2. Every shipped module imports in a fresh interpreter. Import cycles depend on
   which module is entered first, so this is checked one subprocess per module
   rather than in the ambient session.
"""

import functools
import importlib.util
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

import schubmult

# Mirrors the package's own scan in schubmult/__init__.py. _scripts is research
# code that the top-level __getattr__ deliberately tolerates failures from.
_SKIPPED_DIRS = {"__pycache__", "tests", "docs", "scripts", "build", "_scripts"}

# Top-level module names provided by [project.optional-dependencies]. A module
# that fails only because one of these is absent is not a packaging bug, so the
# relevant tests skip instead of failing on installs without the extra.
_OPTIONAL_DEPENDENCIES = frozenset({"matplotlib", "sage"})

_MISSING_MODULE = re.compile(r"No module named '([^']+)'")


def _missing_optional_dependency(message):
    """Return the optional dependency the message blames, if it is genuinely absent."""
    for match in _MISSING_MODULE.finditer(message):
        root = match.group(1).split(".")[0]
        if root in _OPTIONAL_DEPENDENCIES and importlib.util.find_spec(root) is None:
            return root
    return None

# Subpackages that declare a public surface via __all__.
_SUBPACKAGES = [
    "schubmult.combinatorics",
    "schubmult.mult",
    "schubmult.rings",
    "schubmult.rings.combinatorial",
    "schubmult.rings.free_algebra",
    "schubmult.rings.polynomial_algebra",
    "schubmult.rings.schubert",
    "schubmult.symbolic",
]


@functools.lru_cache(maxsize=1)
def _shipped_modules():
    """Every importable module in the installed package, as dotted names."""
    root = Path(schubmult.__file__).resolve().parent
    names = set()
    for py in root.rglob("*.py"):
        parts = py.relative_to(root.parent).with_suffix("").parts
        if any(p in _SKIPPED_DIRS for p in parts):
            continue
        if parts[-1] == "__init__":
            parts = parts[:-1]
        if parts:
            names.add(".".join(parts))
    return sorted(names)


def _subpackage_exports():
    """(package, name) for every name in a subpackage's __all__."""
    import importlib

    pairs = []
    for pkg in _SUBPACKAGES:
        try:
            mod = importlib.import_module(pkg)
        except Exception:  # noqa: BLE001 - reported by the isolation test
            continue
        pairs.extend((pkg, name) for name in getattr(mod, "__all__", []))
    return pairs


def _probe(module):
    result = subprocess.run(
        [sys.executable, "-c", f"import {module}"],
        capture_output=True,
        text=True,
        check=False,
    )
    return module, (result.returncode, result.stderr.strip())


@pytest.fixture(scope="session")
def isolation_results():
    """Import every shipped module in its own interpreter, in parallel."""
    modules = _shipped_modules()
    workers = min(16, (os.cpu_count() or 4) * 2)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return dict(pool.map(_probe, modules))


@pytest.mark.parametrize("module", _shipped_modules())
def test_module_imports_in_isolation(module, isolation_results):
    """A module must import when it is the entry point, not just incidentally.

    Regression guard for the 4.2.0 cycle: schubert_ring imported
    quantum_schubert_ring at load time, so entering at the parabolic module
    raised ImportError on a partially initialized module.
    """
    returncode, stderr = isolation_results[module]
    if returncode != 0:
        missing = _missing_optional_dependency(stderr)
        if missing:
            pytest.skip(f"{module} requires the optional dependency {missing}")
    assert returncode == 0, f"`import {module}` failed in a fresh interpreter:\n{stderr}"


@pytest.mark.parametrize("name", sorted(schubmult._lazy_exports))
def test_declared_top_level_export_resolves(name):
    """Names in _lazy_exports are the documented top-level API.

    __getattr__ turns an underlying ImportError into AttributeError, so a
    broken module surfaces only on attribute access.
    """
    assert getattr(schubmult, name) is not None


@pytest.mark.parametrize(("package", "name"), _subpackage_exports(), ids=lambda v: v.rsplit(".", 1)[-1])
def test_subpackage_export_resolves(package, name):
    import importlib

    assert getattr(importlib.import_module(package), name) is not None


def test_discovered_names_resolve():
    """Names reachable by attribute access must not raise.

    The AST fallback in __getattr__ advertises every public class in every
    module, including modules that cannot be imported. Such a name appears in
    dir(schubmult) and fails only when accessed.
    """
    schubmult._scan_modules()
    broken = {}
    for name, module in sorted(schubmult._module_map.items()):
        try:
            getattr(schubmult, name)
        except Exception as exc:  # noqa: BLE001 - collecting failures to report together
            if _missing_optional_dependency(str(exc)):
                continue
            broken[name] = f"{module}: {exc}"
    assert not broken, "names advertised by dir(schubmult) that fail on access:\n" + "\n".join(f"  {k} <- {v}" for k, v in broken.items())


def test_unknown_attribute_raises_attribute_error():
    """The lazy loader must not mask a typo as something importable."""
    with pytest.raises(AttributeError):
        schubmult.ThisNameDoesNotExist


_NO_SAGE_SCRIPT = r"""
import importlib, importlib.util, sys
# make any `import sage...` fail loudly, whether or not Sage is installed
sys.modules["sage"] = None
import schubmult
for name in schubmult._lazy_exports:
    getattr(schubmult, name)
for module in sys.argv[1:]:
    try:
        importlib.import_module(module)
    except ModuleNotFoundError as exc:
        # other optional extras (matplotlib) may be absent; Sage must never be the reason
        if exc.name is None or exc.name.split(".")[0] in ("sage", "schubmult") or importlib.util.find_spec(exc.name.split(".")[0]) is not None:
            raise
from schubmult import DSx, QPSx, Sx, uncode
assert Sx([3, 1, 2]) * Sx([2, 1, 3])
assert DSx([3, 1, 2]) * DSx([2, 3, 1], "z")
assert QPSx(2, 3)(uncode([2, 3])) * QPSx(2, 3)(uncode([0, 1]))
from schubmult._scripts.schubmult_double import main
main(["schubmult_double", "3", "1", "2", "-", "2", "1", "3", "--no-print"])
print("ok")
"""


def test_package_is_fully_usable_without_sage():
    """`schubmult.sage` is an optional integration: the rest of the package must never import Sage.

    Runs a fresh interpreter with ``sys.modules["sage"] = None`` (so any Sage import raises), imports
    every shipped module outside ``schubmult.sage`` (modules needing another absent optional extra,
    e.g. matplotlib, are tolerated), resolves the whole lazy top-level API, and runs single, double,
    and parabolic quantum products.
    """
    modules = [m for m in _shipped_modules() if not m.startswith("schubmult.sage")]
    result = subprocess.run([sys.executable, "-c", _NO_SAGE_SCRIPT, *modules], capture_output=True, text=True, check=False)
    assert result.returncode == 0 and result.stdout.strip() == "ok", f"schubmult pulled in Sage or failed without it:\n{result.stderr}"
