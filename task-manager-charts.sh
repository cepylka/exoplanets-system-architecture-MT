#!/bin/bash

ls -L1 ./data | while read pickle || [[ -n $pickle ]];
do
    python ./one-system-charts_hz.py "$pickle"
    if [ $? -eq 0 ]; then
        echo -e "Done\n"
    else
        echo "Processing of $pickle failed" >&2
        exit 1
    fi
done
