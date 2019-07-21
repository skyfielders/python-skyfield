#!/bin/bash

set -e
cd "$(dirname ${BASH_SOURCE[0]})"
cd skyfield/data
rm -f Leap_Second.dat deltat.data deltat.preds
python -c 'from skyfield.api import load; load.timescale()' || exit $?
cp Leap_Second.dat deltat.data deltat.preds ../../ci
