#!/bin/bash

set -e
cd "$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"
cd ..

script='import numpy as np
from skyfield.api import load, wgs84
ts = load.timescale()
t = ts.utc(2020, 1, range(365 * 10))
place = wgs84.latlon(0, 0)
v = place.at(t).xyz.km
np.savez_compressed("tmp_before", v=v)'

git restore skyfield/data/iers.npz

echo Before:
python -m skyfield
python -c "$script"
echo

#rm -f finals2000A.all
python builders/build_arrays.py

echo
echo After:
python -m skyfield

script='import numpy as np
with np.load("tmp_before.npz") as t:
    v1 = t["v"]

from skyfield.api import load, tau, wgs84
from skyfield.functions import angle_between
ts = load.timescale()
t = ts.utc(2020, 1, range(365 * 10))
place = wgs84.latlon(0, 0)
v2 = place.at(t).xyz.km
a = angle_between(v1, v2)
a *= 360 * 60 * 60 / tau

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(t.J, a)
ax.grid()
fig.savefig("tmp.png")

print(f"Min {min(a)} Max {max(a)} arcseconds")'

echo
echo Difference:
python -c "$script"
echo '(see tmp.png for a diagram)'
