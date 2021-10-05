#!/usr/bin/env python

from __future__ import print_function

import argparse
import datetime as dt
import sys
from math import atan, degrees, radians, tan

def main(argv):
    parser = argparse.ArgumentParser(description='Build magnitude tests.')
    parser.parse_args(argv)

    # The test data file comes from the paper:
    # https://arxiv.org/pdf/1808.01973.pdf

    with open('Ap_Mag_Input_V3.txt') as f:
        lines = f.read().splitlines()

    print("""
import numpy as np
from numpy import nan
from skyfield import magnitudelib as m
from skyfield.tests.conventions import A""")

    tests = {}  # 'planet': [(args, answer)]

    for line in lines:
        if line.isupper():
            planet = line.lower()
            if planet not in tests:
                tests[planet] = []
        elif line[0] == ' ' and line[1].isdigit():
            fields = line.split()
            if planet in ('mercury', 'earth', 'jupiter'):
                args = fields[4], fields[6], fields[8]
            elif planet == 'venus':
                args = fields[4], fields[6], fields[10]
            elif planet == 'mars':
                args = fields[10], fields[12], fields[16]
            elif planet == 'uranus':
                args = fields[8], fields[10], fields[12], fields[7], fields[5]
            elif planet == 'neptune':
                datestr = fields[0]
                d = dt.datetime.strptime(datestr, '%Y-%b-%d')
                # Formula from Ap_Mag_V3.f90
                year = '%.4f' % (d.year + (d.month - 1) / 12.0 + d.day / 365.0)
                args = fields[8], fields[10], fields[12], year
            else:
                continue
            answer = fields[-1]
            tests[planet].append((args, answer))
        elif line[0] == ' ' and line[1] in 'TF' and planet == 'saturn':
            fields = line.split()

            # The HORIZONS input supplies geodetic latitudes, but the
            # routine needs geocentric latitudes (except that "geo" here
            # is Saturn).
            s = float(fields[6])
            e = float(fields[8])
            e2 = 0.8137e0  # eccentricity squared of Saturn ellipse
            sun_sub_lat = degrees(atan(e2 * tan(radians(s))))
            earth_sub_lat = degrees(atan(e2 * tan(radians(e))))

            args = (
                fields[9], fields[11], fields[15],
                str(sun_sub_lat), str(earth_sub_lat), str(fields[0] == 'T'),
            )
            answer = fields[16] if len(fields) > 16 else 'nan'
            tests[planet].append((args, answer))

    for planet, test_list in tests.items():
        tolerance = (
            # Mars rotation effects are not yet written up.
            '0.1' if planet == 'mars'
            # Other planets should match to high precision.
            else '0.0005'
        )

        if not test_list:
            continue
        print('\ndef test_{}_magnitude_function():'.format(planet))
        for args, answer in test_list:
            joined = ', '.join(str(arg) for arg in args)
            print(f'    mag = m._{planet}_magnitude({joined})')
            if answer == 'nan':
                print(f'    assert np.isnan(mag)')
            else:
                print(f'    assert abs({answer} - mag) < {tolerance}')

        print()
        print('    args = [')
        for arg_vector in zip(*[args for args, answer in test_list]):
            joined = ', '.join(arg_vector)
            print(f'        A[{joined}],')
        print('    ]')
        print(f'    magnitudes = m._{planet}_magnitude(*args)')
        joined = ', '.join(answer for args, answer in test_list)
        print(f'    expected = [{joined}]')
        print(f'    np.allclose(magnitudes, expected, 0, {tolerance},'
              ' equal_nan=True)')

if __name__ == '__main__':
    main(sys.argv[1:])
