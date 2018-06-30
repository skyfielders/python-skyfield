#!/bin/bash

set -e
cd "$(dirname ${BASH_SOURCE[0]})"
make -C skyfield/documentation doctest
