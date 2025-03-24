#!/bin/bash
cat test_cmd_lines.json | xargs -I FILE `FILE` -g > `echo FILE | sed -r 's/--//g' | sed -r 's/ - /T/' | sed -r 's/ /_/g'`
