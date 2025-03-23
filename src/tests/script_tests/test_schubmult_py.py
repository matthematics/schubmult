import sys
def get_json(file: str):
    import os
    import json
    script_dir = os.path.dirname(__file__)
    rel_path = f"data/{file}.json"
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "r") as f:
        return json.load(f)


def assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices):
    from schubmult.schubmult_py import schubmult, schub_coprod
    if coprod:
        coeff_dict = schub_coprod(v_tuple, indices)
    else:
        coeff_dict = schubmult(input_dict, v_tuple)
    for k, v in coeff_dict.items():
        if v != 0:
            assert k in ret_dict
            assert coeff_dict[k] == ret_dict[k]
        else:
            assert v > 0
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert coeff_dict[k] == ret_dict[k]

# def assert_coeff_equal(file_key: str):
#     from ast import literal_eval
#     from schubmult.schubmult_py._script import main
#     from schubmult.perm_lib import uncode
#     input_data = get_json(file_key)
#     coprod = input_data.get("coprod", False)
#     indices = input_data.get("indices", None)
#     ret_dict = main(input_data["cmd_line"] + ["--display-mode", "raw"])
#     ascode = input_data.get("code", False)
#     if ascode:
#         ret_dict = {tuple(uncode(k)): v for k, v in ret_dict.items()}
#     input_dict = {}
#     if not coprod:
#         input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}
#     if input_data.get("direct_check", False):
#         assert input_dict == ret_dict
#     else:        
#         v_tuple = tuple(input_data["v_tuple"])
#         assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices)

# def test_coeff_equal():    
#     json_files = ["schubmult_py_coeff_test1", "schubmult_py_coeff_test2", "schubmult_py_coprod_test1", "really_execute"]
#     for json_file in json_files:
#         assert_coeff_equal(json_file)


def test_coeff_equal_exec(capsys):    
    json_files = ["schubmult_py_coeff_test1", "schubmult_py_coeff_test2", "schubmult_py_coprod_test1", "really_execute"]
    for json_file in json_files:
        from ast import literal_eval
        from schubmult.perm_lib import uncode, permtrim    
        input_data = get_json(json_file)
        from schubmult.schubmult_py._script import main
        # sys.argv = input_data["cmd_line"]
        if "input_dict" in input_data:
            input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}    
        else:
            input_dict = None
        ret_dict = {}
        ascode = input_data.get("code", False)
        coprod = input_data.get("coprod", False)
        indices = input_data.get("indices", None)
        
        main(input_data["cmd_line"])
        lines = capsys.readouterr()
        
        lines = str(lines.out).split("\n")

        if not coprod:        
            for line in lines:            
                try:
                    v, k = line.split("  ")
                except ValueError:
                    print(f"{line=}")
                    continue
                try:
                    v = int(v)
                except ValueError:
                    continue
                ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = int(v)
        else:
            for line in lines:
                line = str(line)
                first_split = ") (" if not ascode else "] ["
                second_split = " (" if not ascode else " ["
                jn = "),(" if not ascode else "],["
                try:
                    vf, s = line.split(first_split)
                    print(f"{(vf, s)}")
                    v, f = vf.split(second_split)
                    print(f"{(v, f)}")
                    k1, k2 = literal_eval(f"({second_split}{f+jn+s})")
                    print(f"{(k1, k2)}")
                    if ascode:
                        k1 = tuple(permtrim(uncode(k1)))
                        k2 = tuple(permtrim(uncode(k2)))
                    k = (k1, k2)
                except ValueError:
                    import traceback
                    print(f"{line=} {first_split=}")
                    print(traceback.format_exc())
                    continue
                v = int(v)
                ret_dict[k] = v

        v_tuple = tuple(input_data["v_tuple"])
        assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices)


if __name__ == "__main__":
    test_coeff_equal_exec()
