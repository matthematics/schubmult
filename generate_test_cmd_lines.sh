#!/bin/bash
find . -wholename '*/script_tests/*.json' | xargs -I FILE python3.10 src/tests/script_tests/_cmd_line_json.py FILE > test_cmd_lines.dat
