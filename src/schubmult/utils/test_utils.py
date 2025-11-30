import symengine


def generate_all(module, filename):
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


def get_json(file: str):
    import json
    import os

    script_dir = os.path.dirname(__file__)
    rel_path = f"../../../tests/script_tests/data/{file}.json"
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path) as f:
        return json.load(f)


def load_json_test_names(this_dir):
    import os

    script_dir = os.path.dirname(__file__)
    rel_path = f"../../../tests/script_tests/data/{this_dir}"
    abs_path = os.path.join(script_dir, rel_path)
    files = os.listdir(abs_path)
    ret = []
    for file in files:
        index = file.rfind(".json")
        filename = file[:index]
        ret += [filename]
    return ret


def print_args(poly):
    def _pr(ag):
        if hasattr(ag, "__sympy__") and not ag.is_Atom:
            return f"({type(ag)},{print_args(ag)})"
        return str(type(ag))

    return "[" + ",".join([_pr(arg) for arg in poly.args]) + "]"


def sympify_args(poly):
    try:
        return symengine.sympify(poly)
    except Exception:
        if poly.is_Mul:
            return symengine.Mul(*[sympify_args(arg) for arg in poly.args])
        if poly.is_Pow:
            return symengine.Pow(*[sympify_args(arg) for arg in poly.args])
        if poly.is_Add:
            return symengine.Add(*[sympify_args(arg) for arg in poly.args])
