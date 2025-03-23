from schubmult.schubmult_py._script import main
from schubmult.schubmult_py import schubmult
from schubmult.perm_lib import uncode
from ast import literal_eval
import os
import sys
import json

def get_json(file: str):
    script_dir = os.path.dirname(__file__)
    rel_path = f"data/{file}.json"
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "r") as f:
        print(f"loading {abs_file_path=}", file = sys.stderr)
        return json.load(f)


def assert_dict_good(v_tuple, input_dict, ret_dict):
    coeff_dict = schubmult(input_dict, v_tuple)
    for k, v in coeff_dict.items():
        if v != 0:
            assert k in ret_dict
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert v > 0
            assert k not in ret_dict
    
    for k in ret_dict.keys():
        assert k in coeff_dict
        assert coeff_dict[k] == ret_dict[k]

def assert_coeff_equal(file_key: str):
    input_data = get_json(file_key)
    ret_dict = main(input_data["cmd_line"] + ["--display-mode", "raw"])
    ascode = input_data.get("code", False)
    if ascode:
        ret_dict = {tuple(uncode(k)): v for k, v in ret_dict.items()}
    input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}
    if input_data.get("direct_check", False):
        assert input_dict == ret_dict
    else:        
        v_tuple = tuple(input_data["v_tuple"])
        assert_dict_good(v_tuple, input_dict, ret_dict)

from io import StringIO
import contextlib

@contextlib.contextmanager
def stdoutIO(stdout=None):
    old = sys.stdout
    if stdout is None:
        stdout = StringIO()
    sys.stdout = stdout
    yield stdout
    sys.stdout = old


def assert_coeff_equal_exec(file_key: str):
    input_data = get_json(file_key)
    sys.argv = input_data["cmd_line"]
    input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}    
    ret_dict = {}
    ascode = input_data.get("code", False)
    with stdoutIO() as s:
        script_dir = os.path.dirname(__file__)
        abs_file_path = os.path.join(script_dir, sys.argv[0])
        print(f"Executing {abs_file_path=}", file=sys.stderr)
        with open(abs_file_path, "r") as script:
            exec(script.read(), globals())
        st = s.getvalue()
        # print(st,file=sys.stderr)
        lines = st.split("\n")

        for line in lines:
            try:
                v, k = line.split("  ")
            except ValueError:
                continue
            v = int(v)
            ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = int(v)
    v_tuple = tuple(input_data["v_tuple"])
    assert_dict_good(v_tuple, input_dict, ret_dict)

        # ret_dict = main(input_data["cmd_line"])
    

#     input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}
#     if input_data.get("direct_check", False):
#         assert input_dict == ret_dict
#     else:        
#         v_tuple = tuple(input_data["v_tuple"])
#         assert_dict_good(v_tuple, input_dict, ret_dict)

def test_coeff_equal():
    json_files = ["schubmult_py_coeff_test1", "schubmult_py_coeff_test2", "schubmult_py_coprod_test1", "really_execute"]
    for json_file in json_files:
        assert_coeff_equal(json_file)
        print(f"Success {json_file}",file=sys.stderr)

def test_coeff_equal_exec():
    json_files = ["schubmult_py_coeff_test1", "schubmult_py_coeff_test2", "really_execute"]
    for json_file in json_files:
        assert_coeff_equal_exec(json_file)
        print(f"Success exec {json_file}",file=sys.stderr)

def print_output(file: str):
    input_data = get_json(file)
    print(main(input_data["cmd_line"]))

if __name__ == "__main__":
    # assert_coeff_equal_exec("really_execute")
    # print_output(sys.argv[1])
    test_coeff_equal()
    test_coeff_equal_exec()