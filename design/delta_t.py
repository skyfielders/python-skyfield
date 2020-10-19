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
from skyfield import functions, timelib
from skyfield.api import load, Loader
from skyfield.data import iers

inf = float('inf')

def main(argv):
    url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt'
    with load.open(url) as f:
        columns = load_table_S15(f)
    i, start_year, end_year, a0, a1, a2, a3 = columns
    report = []

    # Range of years to plot.

    y = np.arange(start_year[0], end_year[-1] + 0.1, 0.01)

    # Skyfield original tables.

    ts = Loader('ci').timescale(builtin=False)
    print('Old shape:', ts.delta_t_table.shape)
    t = ts.J(y)

    T0 = time()
    t.delta_t
    report.append((time() - T0, 's for old interpolation tables'))

    # Skyfield IERS table-driven interpolation.

    ts = load.timescale()
    print('IERS shape:', ts.delta_t_table.shape)
    t = ts.J(y)

    T0 = time()
    old_delta_t = t.delta_t
    report.append((time() - T0, 's for IERS table interpolated delta_t'))

    # New Parabola.

    c1825 = (y - 1825.0) / 100.0
    t0 = time()
    delta_t_parabola = -320.0 + 32.5 * c1825 * c1825
    report.append((time() - t0, 's to compute ∆T with 2018 parabola'))

    # New 2018 splines.

    T0 = time()

    i = np.searchsorted(start_year, y, 'right') - 1
    y0 = start_year[i]
    y1 = end_year[i]
    t = (y - y0) / (y1 - y0)
    delta_t = ((a3[i] * t + a2[i]) * t + a1[i]) * t + a0[i]

    report.append((time() - T0, 's to compute ∆T with 2018 splines'))

    # New 2018 splines combined with IERS data.

    arrays = functions.load_bundled_npy('iers.npz')
    iers_tt, iers_delta_t = arrays['delta_t_recent']
    iers_y = (iers_tt - 1721045.0) / 365.25

    cutoff_index = np.searchsorted(end_year, iers_y[0])
    cutoff_year = end_year[cutoff_index - 1]

    print(cutoff_index, len(end_year))
    print(cutoff_year)
    # print(end_year[:cutoff_index])

    def to_tt(year):
        return year * 365.25 + 1721045.0

    # : i, start_year, end_year, a0, a1, a2, a3 :
    z = np.zeros_like(iers_tt)
    #iers_end_year = 

    cat = np.concatenate
    c = cutoff_index
    start_tt = cat([to_tt(start_year[:c]), to_tt(end_year[-1:]), iers_tt[:-1]])
    end_tt = cat([to_tt(end_year[:c]), iers_tt])
    a0 = cat([a0[:c], [0], iers_delta_t[:-1]])
    a1 = cat([a1[:c], iers_delta_t - a0[-len(iers_delta_t):]])
    a2 = cat([a2[:c], z])
    a3 = cat([a3[:c], z])

    # Try using combined splines.

    T0 = time()

    tt = to_tt(y)
    i = np.searchsorted(start_tt, tt, 'right') - 1
    tt0 = start_tt[i]
    tt1 = end_tt[i]
    t = (tt - tt0) / (tt1 - tt0)
    experimental_delta_t = ((a3[i] * t + a2[i]) * t + a1[i]) * t + a0[i]

    report.append((time() - T0, 's to compute ∆T with combined spline table'))

    # The plot.

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    ax = ax1

    ax.plot(y, delta_t_parabola, label='2018 parabola')
    ax.plot(y, delta_t, label='2018 splines ∆T')
    ax.plot(y, old_delta_t, label='Old Skyfield ∆T')

    #ax.set_xlim(-720, -718)
    #ax.set_ylim(-400, 200)
    ax.set_xlim(1550, 2020)
    ax.set_ylim(-400, 200)

    ax.legend()
    ax.grid()

    ax = ax2

    i = (y >= 1973)
    ax.plot(y[i], delta_t[i] - old_delta_t[i])
    ax.grid()
    ax.set(ylabel='2018 splines - IERS')

    ax = ax3

    end_year = 1973.1
    i = (y >= 1972.9) & (y <= end_year)
    ax.plot(y[i], delta_t[i])

    i = (iers_y <= end_year)
    ax.plot(iers_y[i], iers_delta_t[i], '.')
    ax.grid()
    #ax.set(ylabel='2018 splines - IERS')

    fig.savefig('tmp.png')

    for args in report:
        print(*args)

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

