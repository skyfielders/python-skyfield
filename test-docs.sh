#!/bin/bash

set -e
cd "$(dirname ${BASH_SOURCE[0]})"

# Copy in files that would otherwise need to download.
cp ci/* skyfield/documentation

function cleanup {
    cd skyfield/documentation
    rm -f $(cd ../../ci; ls)
}
trap cleanup EXIT

make -C skyfield/documentation doctest

# Remove the files we copied in.
FILES="$(cd ci; ls *.*)"
(cd skyfield/documentation; rm $FILES)
