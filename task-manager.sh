#!/bin/bash

cat systems-min2planets-withMassSinAndRadius.txt | while read star || [[ -n $star ]];
do
   python ./one-system_rewrite_star.py "$star"
   if [ $? -eq 0 ]; then
        echo -e "Done\n"
    else
        echo "Processing of $star failed" >&2
        exit 1
    fi
done
