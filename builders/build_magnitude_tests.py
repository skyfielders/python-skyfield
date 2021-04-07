#!/usr/bin/env python

from __future__ import print_function

import argparse
import sys

def main(argv):
    parser = argparse.ArgumentParser(description='Build magnitude tests.')
    parser.parse_args(argv)

    # The test data file comes from the paper:
    # https://arxiv.org/pdf/1808.01973.pdf

    with open('Ap_Mag_Input_V3.txt') as f:
        lines = f.read().splitlines()

    print("""
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
            elif planet == 'uranus':
                args = fields[8], fields[10], fields[12], fields[7], fields[5]
            else:
                continue
            answer = fields[-1]
            tests[planet].append((args, answer))

    # from pprint import pprint
    # pprint(tests)

    for planet, test_list in sorted(tests.items()):
        if not test_list:
            continue
        print('\ndef test_{}_magnitude_function():'.format(planet))
        for args, answer in test_list:
            joined = ', '.join(str(arg) for arg in args)
            print(f'    mag = m._{planet}_magnitude({joined})')
            print(f'    assert abs({answer} - mag) < 0.0005')

        if planet != 'venus':
            continue

        print()
        print('    args = [')
        for arg_vector in zip(*[args for args, answer in test_list]):
            joined = ', '.join(arg_vector)
            print(f'        A[{joined}],')
        print('    ]')
        print(f'    magnitudes = m._{planet}_magnitude(*args)')
        joined = ', '.join(answer for args, answer in test_list)
        print(f'    expected = [{joined}]')
        print('    assert all(magnitudes - expected < 0.0005)')

if __name__ == '__main__':
    main(sys.argv[1:])
