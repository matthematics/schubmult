"""Helpers for the test suite: locating JSON test data and inspecting SymPy/SymEngine expression trees."""

import symengine


def generate_all(module, filename):
    """Print an import block and ``__all__`` list for the public names defined in ``filename`` (dev helper)."""
    D = dir(module)
    print(f"{D=}")
    file_data = ""
    with open(filename) as f:
        file_data = str(f.read())
    print(f"{file_data}")
    real_d = [d for d in D if (file_data.find(f"def {d}") != -1 or file_data.find(f"class {d}") != -1) and d[0] != "_"]
    print("from bob import (")
    print("    ", end="")
    print(",\n    ".join(real_d))
    print(",")
    print(")")

    print("__all__ =")
    print("[")
    print("    ", end="")
    print("',\n    '".join(real_d))
    print(",")
    print("]")


def _data_dir():
    """tests/scripts/data, located from the calling test module (works for non-editable installs too)."""
    import inspect
    import os

    for frame in inspect.stack()[2:]:
        d = os.path.dirname(os.path.abspath(frame.filename))
        while d and d != os.path.dirname(d):
            cand = os.path.join(d, "tests", "scripts", "data")
            if os.path.isdir(cand):
                return cand
            if os.path.basename(d) == "scripts" and os.path.isdir(os.path.join(d, "data")):
                return os.path.join(d, "data")
            d = os.path.dirname(d)
    raise FileNotFoundError("tests/scripts/data not found relative to the calling test module")


def get_json(file: str):
    """Load ``<file>.json`` from the test data directory."""
    import json
    import os

    with open(os.path.join(_data_dir(), f"{file}.json")) as f:
        return json.load(f)


def load_json_test_names(this_dir):
    """Names (without ``.json``) of all test case files in the data subdirectory ``this_dir``."""
    import os

    files = os.listdir(os.path.join(_data_dir(), this_dir))
    ret = []
    for file in files:
        index = file.rfind(".json")
        filename = file[:index]
        ret += [filename]
    return ret


def print_args(poly):
    """Nested string of the argument types of an expression tree (for debugging printing issues)."""
    def _pr(ag):
        if hasattr(ag, "__sympy__") and not ag.is_Atom:
            return f"({type(ag)},{print_args(ag)})"
        return str(type(ag))

    return "[" + ",".join([_pr(arg) for arg in poly.args]) + "]"


def sympify_args(poly):
    """Convert a SymPy expression to SymEngine, recursing into ``Mul``/``Pow``/``Add`` when direct
    conversion fails.
    """
    try:
        return symengine.sympify(poly)
    except Exception:
        if poly.is_Mul:
            return symengine.Mul(*[sympify_args(arg) for arg in poly.args])
        if poly.is_Pow:
            return symengine.Pow(*[sympify_args(arg) for arg in poly.args])
        if poly.is_Add:
            return symengine.Add(*[sympify_args(arg) for arg in poly.args])
