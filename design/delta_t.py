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

import matplotlib.pyplot as plt
import numpy as np
from skyfield.api import load
from skyfield.data import iers

def main(argv):
    # parser = argparse.ArgumentParser(description=put description here)
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    # args = parser.parse_args(argv)
    # print(args.accumulate(args.integers))

    draw_plot_comparing_USNO_and_IERS_data()

RE = re.compile(b'^(..)(..)(..)' + b'.' * 52 + b'(.\d........)', re.M)

def draw_plot_comparing_USNO_and_IERS_data():
    f = load.open('finals.all')
    year, month, day, dut1 = iers.parse_dut1_from_finals_all(f)

    ts = load.timescale()
    t = ts.utc(year, month, day)
    y_new = dut1
    y_old = t.dut1

    fig, ax = plt.subplots()
    ax.set(xlabel='Year', title='UT1 minus UTC')
    ax.plot(t.J, y_old, '-', label='deltat.dat')
    ax.plot(t.J, y_new, ',', label='finals.all')
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

