#!/bin/bash

set -e
cd "$(dirname ${BASH_SOURCE[0]})"

# Copy in files that would otherwise need to download.
ln -f ci/*.* skyfield/documentation

function cleanup {
    cd skyfield/documentation
    rm -f $(cd ../../ci; ls *.*)
}
trap cleanup EXIT

make -C skyfield/documentation doctest
