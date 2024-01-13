#!/bin/bash

set -e
cd "$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"

echo Before:
python -m skyfield
echo

cd ..
rm -f finals2000A.all
python builders/build_arrays.py

echo
echo After:
python -m skyfield
