#!/usr/bin/env python
"""Print the code for the skyfield/data/spice.py file.

"""
import re
import os
import sys
from pprint import pformat
from textwrap import dedent

from skyfield.constants import ASEC2RAD
from skyfield.functions import rot_x, rot_y, rot_z

axes = {'1': rot_x, '2': rot_y, '3': rot_z}

template = """\
# Machine generated - see build_spice.py

from numpy import array

inertial_frames = [
    %s,
    ]

inertial_frames = dict((key, array(value)) for key, value in inertial_frames)
"""

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: build_spice_rotations.py /path/to/spice/toolkit/",
              file=sys.stderr)
        sys.exit(2)
    path = os.path.join(sys.argv[1], 'src/spicelib/chgirf.f')
    with open(path) as f:
        fortran = f.read()
    fortran = re.sub(r'\n +\.', '', fortran, flags=re.S)
    names = re.findall(r"DATA FRAMES[^']*'([^']+)", fortran)
    bases = re.findall(r"DATA BASES[^']*'([^']+)", fortran)
    rotations = re.findall(r"DATA DEFS[^']*'([^']+)", fortran)
    assert len(names) == len(bases) == len(rotations) == 21
    bases = {name: base for name, base in zip(names, bases)}
    frames = {}
    framelist = []
    for name, rotation in zip(names, rotations):
        fields = rotation.replace(',', ' ').split()
        M = None
        for i in reversed(range(0, len(fields), 2)):
            theta, axis = fields[i:i+2]
            theta = float(theta.replace('D', 'E')) * ASEC2RAD
            f = axes[axis]
            R = f(-theta)
            M = R if (M is None) else R * M
        base = bases[name]
        while base != 'J2000':
            M = frames[base] * M
            base = bases[base]
        lists = list(list(row) for row in M)
        frames[name] = lists
        framelist.append((name, lists))

    # The pformat() function produces terrible output if its own
    # `indent` keyword is set, so we rig up our own indentation, and
    # also provide a trailing comma:

    print(template
          % dedent(
              ' ' + pformat(framelist).strip('[]')
          ).replace('\n', '\n    '), end='')
