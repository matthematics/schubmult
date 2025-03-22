from schubmult.schubmult_py._script import main, schubmult
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


def assert_coeff_equal(file_key: str):
    input_data = get_json(file_key)
    ret_dict = main(input_data["cmd_line"])
    input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}
    if input_data.get("direct_check", False):
        assert input_dict == ret_dict
    else:        
        v_tuple = tuple(input_data["v_tuple"])
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

def test_coeff_equal():
    json_files = ["schubmult_py_coeff_test1", "schubmult_py_coeff_test2", "schubmult_py_coprod_test1"]
    for json_file in json_files:
        assert_coeff_equal(json_file)
        print("Success",file=sys.stderr)

def print_output(file: str):
    input_data = get_json(file)
    print(main(input_data["cmd_line"]))

if __name__ == "__main__":
    test_coeff_equal()
    # print_output(sys.argv[1])