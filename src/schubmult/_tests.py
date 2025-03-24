
def get_json(file: str):
    import json
    import os

    script_dir = os.path.dirname(__file__)
    rel_path = f"../tests/script_tests/data/{file}.json"
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path) as f:
        return json.load(f)


def load_json_test_names(this_dir):
    import os
    script_dir = os.path.dirname(__file__)
    rel_path = f"../tests/script_tests/data/{this_dir}"
    abs_path = os.path.join(script_dir, rel_path)
    files = os.listdir(abs_path)
    ret = []
    for file in files:
        index = file.rfind(".json")
        filename = file[:index]
        ret += [filename]
    return ret
