#!/bin/bash

# help
if [ $# -lt 3 ]; then
    echo "Usage: $0 file_prefix endnum ext"
    echo "Ex: $0 chr6_haps_ 9 .v.maxl"
    echo "This will check for chr6_haps_0.v.maxl to chr6_haps_9.v.maxl"
    exit 1
fi

# args
prefix=$1
end=$2
suffix=$3

# start at 0
start=0

for i in $(seq $start $end); do
    filename="${prefix}${i}${suffix}"
    if [ ! -f "$filename" ]; then
        echo "Missing: $filename"
    fi
done