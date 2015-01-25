"""Rebuild the test data that compares Skyfield to the NOVAS library."""

from __future__ import print_function

import sys
from datetime import datetime
from itertools import product
from textwrap import dedent

banner = '*' * 79
one_second = 1.0 / 24.0 / 60.0 / 60.0
planets = [('mercury', 1), ('venus', 2), ('mars', 4), ('jupiter', 5),
           ('saturn', 6), ('uranus', 7), ('neptune', 8), ('pluto', 9),
           ('sun', 10), ('moon', 11)]

def main(in_path):
    output({}, """\
        'Auto-generated accuracy tests vs HORIZONS (build_horizons_tests.py).'

        from numpy import max
        from skyfield import api
        from skyfield.constants import AU_M

        one_second = 1.0 / 24.0 / 60.0 / 60.0
        arcsecond = 1.0 / 60.0 / 60.0
        ra_arcsecond = 24.0 / 360.0 / 60.0 / 60.0
        meter = 1.0 / AU_M

        def compare(value, expected_value, epsilon):
            if hasattr(value, 'shape') or hasattr(expected_value, 'shape'):
                assert max(abs(value - expected_value)) <= epsilon
            else:
                assert abs(value - expected_value) <= epsilon

        """)

    lines = read_lines(open(in_path))
    bannered = False
    i = 1

    for line in lines:
        #print(repr(line))
        if bannered:
            field_names = line
            bannered = False
        if line.startswith('Target body name: '):
            name = line.split()[3].lower()
        elif line == banner:
            bannered = True
        elif line == '$$SOE':
            while True:
                line = next(lines)
                if line == '$$EOE':
                    break
                generate_test(i, name, field_names, line)
                i += 1


def generate_test(i, name, field_names, line):
    date = get_date(field_names, line)
    hlon = get_float(field_names, line, 'hEcl-Lon')
    hlat = get_float(field_names, line, 'hEcl-Lat')
    output(locals(), """\
        def test_{name}{i}():
            astrometric = api.sun(utc={date}).observe(api.{name})
            hlat, hlon, d = astrometric.ecliptic_latlon()
            compare(hlat.degrees, {hlat}, 0.001)
            compare(hlon.degrees, {hlon}, 0.001)
        """)


def get_date(field_names, line):
    t = get_text(field_names, line, 'Date__(UT)__HR:MN')
    d = datetime.strptime(t, '%Y-%b-%d %H:%M')
    return (d.year, d.month, d.day, d.hour, d.minute)


def get_float(field_names, line, field_name):
    return float(get_text(field_names, line, field_name))


def get_text(field_names, line, field_name):
    i = field_names.index(field_name)
    return line[i:i+len(field_name)].strip()


def read_lines(file):
    for line in file:
        yield line.rstrip()


def output(dictionary, template):
    print(dedent(template).format(**dictionary).strip('\n'))
    print()


if __name__ == '__main__':
    main(sys.argv[1])
