#!/bin/bash

set -e
cd "$(dirname ${BASH_SOURCE[0]})"

# Copy in files that would otherwise need to download.
ln -f ci/*.* documentation/

function cleanup {
    cd documentation
    rm -f $(cd ../ci && ls *.*)
}
trap cleanup EXIT

make -C documentation doctest
