def get_json(file: str):
    import os
    import json

    script_dir = os.path.dirname(__file__)
    rel_path = f"../tests/{file}.json"
    abs_file_path = os.path.join(script_dir, rel_path)
    with open(abs_file_path, "r") as f:
        return json.load(f)