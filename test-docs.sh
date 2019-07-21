#!/bin/bash

set -e
cd "$(dirname ${BASH_SOURCE[0]})"
cp ci/* skyfield/documentation

function cleanup {
    cd skyfield/documentation
    rm -f $(cd ../../ci; ls)
}
trap cleanup EXIT

make -C skyfield/documentation doctest
