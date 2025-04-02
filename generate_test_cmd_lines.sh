#!/bin/bash
find . -wholename './src/tests/script_tests/data/*/*.json' | xargs -I FILE python3.10 src/tests/script_tests/_cmd_line_json.py FILE > test_cmd_lines.dat
