"""
Choosing a source for ∆T
========================

1. `ftp://ftp.iers.org/products/eop/rapid/standard/`_
2. `ftp://ftp.iers.org/products/eop/rapid/standard/csv/`_

* The CSV files are slightly larger, both compressed and uncompressed.
"""

#!/usr/bin/env python3

import argparse
import re
import sys
from time import time

import matplotlib.pyplot as plt
import numpy as np
from skyfield import timelib
from skyfield.api import load
from skyfield.data import iers

inf = float('inf')

def main(argv):
    # parser = argparse.ArgumentParser(description=put description here)
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    # args = parser.parse_args(argv)
    # print(args.accumulate(args.integers))

    #draw_plot_comparing_USNO_and_IERS_data()
    reduce_IERS_data()

RE = re.compile(b'^(..)(..)(..)' + b'.' * 52 + b'(.\d........)', re.M)

def reduce_IERS_data():
    f = load.open('finals.all')
    mjd, dut1 = iers.parse_dut1_from_finals_all(f)
    jd = mjd + 2400000.5
    # print(dut1)

    leap_second_indices = np.concatenate([[False], np.diff(dut1) > 0.9])
    delta_t2 = np.cumsum(leap_second_indices) - dut1 + 32.184 + 12.0

    delta_t_recent = np.array([jd, delta_t2])
    #print(delta_t_recent)

    leap_dates = 2400000.5 + np.concatenate([
        [-inf], [41317.0, 41499.0, 41683.0], mjd[leap_second_indices], [inf],
    ])
    #print(leap_dates)
    leap_offsets = np.concatenate([
        [10.0, 10.0], np.arange(10, len(leap_dates) + 8),
    ])
    print(leap_offsets)

    ts = load.timescale()
    ts
    #print(ts.delta_t_table)

    jd = mjd + 2400000.5
    t = ts.tt_jd(jd)

    #print(delta_t2 - t.delta_t)

    fig, ax = plt.subplots()
    ax.plot(t.J, t.delta_t)
    ax.plot(t.J, delta_t2)
    ax.grid()
    fig.savefig('tmp.png')

def draw_plot_comparing_USNO_and_IERS_data():
    f = load.open('finals.all')
    year, month, day, dut1 = iers.parse_dut1_from_finals_all(f)

    # auto-detect leap seconds
    leap_seconds = np.diff(dut1) > 0.9
    #x = np.append(leap_seconds > 0.9, False)
    #x = np.prepend(leap_seconds > 0.9, False)
    #x = np.concatenate([leap_seconds > 0.9, [False]]).astype(bool)
    x = np.concatenate([[False], leap_seconds, ])
    print(x)
    print(x, sum(x))
    print(year[x])
    print(month[x])
    print(day[x])

    delta_t = dut1 - np.cumsum(x)

    ts = load.timescale()
    t = ts.utc(year, month, day)
    y_new = dut1
    y_old = t.dut1



    # and: figure out error cost of interpolating weeks or months
    all_points = np.arange(len(delta_t))
    samples = all_points[::600] + 0.001
    print(f'{len(samples)} samples')

    for skip in range(1, 40):
        index = np.arange(0, len(delta_t), skip)
        excerpt = delta_t[index]
        interpolated = np.interp(all_points, index, excerpt)
        difference_seconds = abs(interpolated - delta_t).max()
        difference_arcseconds = difference_seconds / 24.0 * 360.0
        storage = 4 * len(excerpt)

        t0 = time()
        np.interp(samples, index, excerpt)
        duration = time() - t0

        fmt = '{:2} days  {:10.6f} s   {:10.6f} arcseconds  {} bytes  {:.6f} s'
        print(fmt.format(
            skip, difference_seconds, difference_arcseconds, storage, duration,
        ))
    #print(dut1)

    fig, ax = plt.subplots()
    ax.set(xlabel='Year', title='UT1 minus UTC')
    # ax.plot(t.J, y_old, '-', label='deltat.dat')
    # ax.plot(t.J, y_new, ',', label='finals.all')
    ax.plot(delta_t)
    ax.plot(np.diff(delta_t))
    ax.grid()
    plt.legend()
    fig.savefig('tmp.png')

    print(len(y_new))
    np.savez_compressed(
        'test.npz',
        time=t.J,
        ut1_minus_utc=y_new,
    )

"""
Example line from finals.all, for crafting our RE:

73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916  P    44.969     .500     2.839     .300   .143000   .137000   .8075000      .000      .000  \n\\
"""

if __name__ == '__main__':
    main(sys.argv[1:])

