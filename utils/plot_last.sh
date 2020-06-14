#!/bin/bash

file=`(ls -tc ./data | head -n1)`
echo $file
python3 ./utils/plot.py $file
