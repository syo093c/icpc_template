#! /bin/bash -eu

filename=$(basename "$1" | cut -d. -f1)
clang++ $1 -o $filename && ./$filename
