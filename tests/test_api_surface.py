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
