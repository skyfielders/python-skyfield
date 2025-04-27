#!/usr/bin/env python3

import argparse
import sys
from math import pi
from skyfield import almanac, api
from skyfield.api import N, W, load, wgs84

def main(argv):
    parser = argparse.ArgumentParser(description='Print conjunctions')
    parser.add_argument('year', type=int, help='year')
    parser.add_argument('target1', help='name of first target')
    parser.add_argument('target2', help='name of second target')
    args = parser.parse_args(argv)

    ts = load.timescale()
    t0 = ts.utc(args.year, 0, 0)
    t1 = ts.utc(args.year + 1, 0, 0)

    eph = load('de421.bsp')
    earth = eph['earth']
    target1 = eph[args.target1]
    target2 = eph[args.target2]

    def f(t):
        e = earth.at(t)
        _, lon1, _ = e.observe(target1).ecliptic_latlon()
        _, lon2, _ = e.observe(target2).ecliptic_latlon()
        return (lon1.radians - lon2.radians) // pi % 2

    f.step_days = 14
    t, y = almanac.find_discrete(t0, t1, f)
    for ti, yi in zip(t, y):
        e = earth.at(ti)
        p1 = e.observe(target1)
        p2 = e.observe(target2)
        print('{}  {}  {:6.2f}'.format(
            ti.utc_strftime(), yi, p1.separation_from(p2).degrees,
        ))

if __name__ == '__main__':
    main(sys.argv[1:])
