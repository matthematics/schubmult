import json
import os
import sys

with open(sys.argv[1]) as f:
    d = json.load(f)
    dirname = os.path.dirname(sys.argv[1])
    index = dirname.rfind("/")
    filename = dirname[index + 1 :]
    d["cmd_line"][0] = filename
    cmd = " ".join(d["cmd_line"])
    new_name = f"{dirname}/"
    new_name = new_name + cmd.replace("--", "").replace(" - ", "T").replace(" ", "_")
    # print(f"mv {sys.argv[1]} )
    os.rename(sys.argv[1], f"{new_name}.json")
    print(cmd)
