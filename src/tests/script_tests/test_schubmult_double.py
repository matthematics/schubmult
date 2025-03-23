import pytest
from ast import literal_eval
import re


from tests._tests import get_json




def assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same=True):
    print(f"{input_dict=}")
    from schubmult.schubmult_double import schubmult, schub_coprod
    from symengine import sympify, symarray, expand
    var_a = symarray("y", 100)
    var_b = symarray("z", 100)

    if coprod:
        coeff_dict = schub_coprod(v_tuple, indices, var2=var_a, var3=var_a if same else var_b)
    else:
        coeff_dict = schubmult(input_dict, v_tuple, var2=var_a, var3=var_a if same else var_b)

    # print(f"{coeff_dict=}")
    # print(f"{ret_dict=}")

    for k, v in coeff_dict.items():
        #print(f" bef {k=} {v=}")
        #v = expand(sympify(v))
        print(f" {k=} {v=} {ret_dict.get(k, 0)=}")
        assert ((k not in ret_dict) or (expand(v) == expand(ret_dict[k])))        
    for k in ret_dict.keys():
        assert expand(coeff_dict[k]) == expand(ret_dict[k])

def parse_ret(lines, ascode, coprod):
    from schubmult.perm_lib import uncode, permtrim
    from schubmult.schubmult_double._script import main
    from symengine import sympify, symarray, expand
    ret_dict = {}
    if not coprod:
        for line in lines:
            try:
                k, v = line.lstrip().split("  ")
            except Exception as e:                
                print(f"BOOB {coprod=} {line=} {e=}")
                # if len(line)>0:
                #     raise
                continue
            try:
                v = sympify(v)
            except Exception as e:
                print(f"BOOB {v=} {e=}")
                continue
            ret_dict[(tuple(literal_eval(k)) if not ascode else tuple(uncode(literal_eval(k))))] = v
    else:
        nf = 0
        for line in lines:
            line = str(line)
            charo = "[" if ascode else "("
            charc = "]" if ascode else ")"
            first_split = "\\) +\\(" if not ascode else "\\] +\\["
            second_split = ")  " if not ascode else "]  "
            jn = ")," if not ascode else "],"
            try:
                s, vf = re.split(first_split, line)
                print(f"BABO {s=} {vf=}")
                f, v = vf.split(second_split)  # re.split(second_split, vf)
                print(f"BABO {f=} {v=}")
                evlaf = f"({charo}{f + jn + s}{charc})"
                k1, k2 = literal_eval(evlaf)
                print(f"BABO {k1=} {k2=}")
                print(f"BABO {line=}")
                if ascode:
                    k1 = tuple(permtrim(uncode(k1)))
                    k2 = tuple(permtrim(uncode(k2)))
                k = (k2, k1)
            except Exception as e:
                print(f"BOOB foop {line=} {e=}")
                # if len(line)>0:
                #     raise
                # nf += 1
                continue
            try:
                v = sympify(v)
            except Exception as e:
                print(f"BOOB {v=} {e=}")
                continue
            ret_dict[k] = v
            print(f"BABO {k=} {ret_dict[k]=}{v=}")
    return ret_dict

# json_files_data = [
#     "schubmult_double_cf1",
#     "schubmult_double_cf2",
#     "schubmult_double_coprod_mixed",
#     "schubmult_double_coprod_mixed_positive",
# ]


# @pytest.mark.parametrize("json_file", json_files_data)
# def test_coeff_equal_exec(capsys, json_file):
#     import re    
#     from schubmult.perm_lib import uncode, permtrim

#     input_data = get_json(json_file)
#     from schubmult.schubmult_double._script import main
#     from symengine import sympify, symarray, expand

#     # sys.argv = input_data["cmd_line"]
#     ascode = input_data.get("code", False)
#     print(f"BOOB {json_file=} {input_data=}") 
#     input_dict = {
#         literal_eval(k) if not ascode else tuple(permtrim(uncode(literal_eval(k)))): v
#         for k, v in input_data.get("input_dict",{}).items()
#     }    
#     ret_dict = {}
#     coprod = input_data.get("coprod", False)
#     indices = input_data.get("indices", None)
#     same = not input_data.get("mixed", False)
#     main(input_data["cmd_line"])
#     lines = capsys.readouterr()
#     print(f"BOOB {json_file=} {lines=} forple")
#     lines = str(lines.out).split("\n")

#     ret_dict = parse_ret(lines, ascode, coprod)
#     v_tuple = tuple(input_data["v_tuple"]) if not ascode else tuple(uncode(input_data["v_tuple"]))
#     print("BOOB ASS")
#     assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same)
#     print(f"BOOB {json_file=} yay")


base_dir = "script_tests/data/schubmult_double"

json_files_data_args = [f"{base_dir}/test_gen_double", f"{base_dir}/test_gen_double_coprod"]


@pytest.mark.parametrize("json_file", json_files_data_args)
def test_with_same_args_exec(capsys, json_file):
    import re
    from ast import literal_eval
    from schubmult.perm_lib import uncode, permtrim

    args = get_json(json_file)
    print(f"{json_file=} {args=} input_data")
    from schubmult.schubmult_double._script import main
    from symengine import sympify, symarray, expand

    mult = args["mult"]
    mulstring = args["mulstring"]

    perms = args["perms"]

    ascode = args["ascode"]
    coprod = args["coprod"]
    same = args["same"]
    msg = args["msg"]
    down = args["down"]
    display_positive = args["display_positive"]
    pr = args["pr"]

    print(f"{args=}")
    print(f"{args['cmd_line']=}")
    main(args["cmd_line"])
    lines = capsys.readouterr()
    lines = str(lines.out).split("\n")

    ret_dict = parse_ret(lines, ascode, coprod)
    v_tuple = (
        (tuple(perms[1]) if not ascode else tuple(uncode(perms[1])))
        if not coprod
        else (tuple(perms[0]) if not ascode else tuple(uncode(perms[0])))
    )
    input_dict = (
        {tuple(permtrim(perms[0])) if not ascode else tuple(permtrim(uncode(perms[0]))): 1}     
    )
    indices = tuple(perms[1])
    print(f"{v_tuple=} {input_dict=} {indices=}")
    print("BOOB ASS")
    assert_dict_good(v_tuple, input_dict, ret_dict, coprod, indices, same)


if __name__ == "__main__":
    # test_coeff_equal_exec()
    test_with_same_args_exec()
