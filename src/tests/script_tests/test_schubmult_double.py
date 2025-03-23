from symengine import sympify, symarray, expand

var_a = symarray("y", 100)
var_b = symarray("z", 100)

def get_json(file: str):
    import os
    import json
    script_dir = os.path.dirname(__file__)
    rel_path = f"data/{file}.json"
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "r") as f:
        return json.load(f)


def assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same=True):
    from schubmult.schubmult_double import schubmult, schub_coprod
    if coprod:
        coeff_dict = schub_coprod(v_tuple, indices, var2=var_a, var3=var_a if same else var_b)
    else:
        coeff_dict = schubmult(input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b)
    
    for k, v in coeff_dict.items():
        v = expand(v)
        if v != 0:
            assert k in ret_dict
            assert v == expand(ret_dict[k])
        else:            
            assert k not in ret_dict
    for k in ret_dict.keys():
        assert expand(coeff_dict[k]) == expand(ret_dict[k])

def test_coeff_equal_exec(capsys):    
    json_files = ["schubmult_double_cf1", 
                  "schubmult_double_cf2"]
    for json_file in json_files:
        from ast import literal_eval
        from schubmult.perm_lib import uncode, permtrim    
        input_data = get_json(json_file)
        from schubmult.schubmult_double._script import main
        # sys.argv = input_data["cmd_line"]
        if "input_dict" in input_data:
            input_dict = {literal_eval(k): v for k, v in input_data["input_dict"].items()}    
        else:
            input_dict = None
        ret_dict = {}
        ascode = input_data.get("code", False)
        coprod = input_data.get("coprod", False)
        indices = input_data.get("indices", None)
        same = not input_data.get("mixed", False)
        main(input_data["cmd_line"])
        lines = capsys.readouterr()
        lines = str(lines.out).split("\n")

        if not coprod:        
            for line in lines:            
                try:
                    k, v = line.lstrip().split("  ")
                except ValueError:
                    continue
                try:
                    v = sympify(v)
                except Exception:
                    raise
                ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = v
        else:
            for line in lines:
                line = str(line)
                first_split = ") (" if not ascode else "] ["
                second_split = " (" if not ascode else " ["
                jn = "),(" if not ascode else "],["
                try:
                    vf, s = line.split(first_split)
                    v, f = vf.split(second_split)
                    k1, k2 = literal_eval(f"({second_split}{f+jn+s})")
                    if ascode:
                        k1 = tuple(permtrim(uncode(k1)))
                        k2 = tuple(permtrim(uncode(k2)))
                    k = (k1, k2)
                except ValueError:
                    continue
                v = sympify(v)
                ret_dict[k] = v

            v_tuple = tuple(input_data["v_tuple"])
            assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same)


if __name__ == "__main__":
    test_coeff_equal_exec()
