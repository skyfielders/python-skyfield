#!/usr/bin/env python3

"""
Choosing a source for ∆T
========================

1. `ftp://ftp.iers.org/products/eop/rapid/standard/`_
2. `ftp://ftp.iers.org/products/eop/rapid/standard/csv/`_

* The CSV files are slightly larger, both compressed and uncompressed.
* For the position of an object at the celestial equator,
  1 second of time = 15 arcseconds.
  1 millisecond of time = 15 mas.

"""
import argparse
import csv
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
    url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt'
    with load.open(url) as f:
        columns = load_table_S15(f)
    i, start_year, end_year, a0, a1, a2, a3 = columns

    y = np.arange(start_year[0], end_year[-1] + 0.1, 0.1)
    i = np.searchsorted(start_year, y, 'right') - 1
    y0 = start_year[i]
    y1 = end_year[i]
    t = (y - y0) / (y1 - y0)

    from time import time
    t0 = time()
    delta_t = ((a3[i] * t + a2[i]) * t + a1[i]) * t + a0[i]
    print(time() - t0)

    fig, ax = plt.subplots()
    ax.plot(y, delta_t)
    ax.set_xlim(1550, 2020)
    ax.set_ylim(-60, 200)
    ax.grid()
    fig.savefig('tmp.png')

    return

    # parser = argparse.ArgumentParser(description=put description here)
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    # args = parser.parse_args(argv)
    # print(args.accumulate(args.integers))

    #draw_plot_comparing_USNO_and_IERS_data()

    f = load.open('finals2000A.all')
    mjd_utc, dut1 = iers.parse_dut1_from_finals_all(f)
    delta_t_recent, leap_dates, leap_offsets = (
        iers.build_timescale_arrays(mjd_utc, dut1)
    )

    ts = load.timescale()

    print('TT 2050:', ts.J(2050).tt)

    fig, ax = plt.subplots()

    # jd = mjd_utc + 2400000.5
    # t = ts.tt_jd(jd)
    # ax.plot(t.J, t.delta_t)

    t = ts.J(range(1600, 2150))
    ax.plot(t.J, t.delta_t, label='Skyfield')

    t = ts.J(range(2150, 2190))
    ms2004 = timelib.delta_t_formula_morrison_and_stephenson_2004(t.tt)
    ax.plot(t.J, ms2004, label='M&S 2004 long-term parabola')

    y = np.arange(2005, 2050)
    t = y - 2000
    ax.plot(y, 62.92 + 0.32217 * t + 0.005589 * t * t,
            label='M&S 2004 polynomial 2005-2050')

    y = np.arange(2050, 2150)
    ax.plot(y, -20 + 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y),
            label='M&S 2004 polynomial 2050-2150')
    # y = np.array([2050, 2149])
    # ax.plot(y, -20 + 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y))

    ax.set(xlabel='Year', ylabel='∆T')
    ax.grid()
    ax.legend()
    fig.savefig('tmp.png')

    is_2050_to_2150_polynomial_worth_it()

def load_table_S15(f):
    # http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt
    content = f.read()
    banner = b'- ' * 36 + b'-\n'
    sections = content.split(banner)
    table = np.loadtxt(sections[2].splitlines())
    return table.T

def is_2050_to_2150_polynomial_worth_it():
    # How different is the polynomial at:
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    # from simply doing a linear interpolation from the endpoints?
    y = np.arange(2050, 2150)
    delta_t = -20 + 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y)
    linear_approximation = np.linspace(delta_t[0], delta_t[-1], len(delta_t))
    print('Maximum difference:', max(abs(delta_t - linear_approximation)), 's')

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

