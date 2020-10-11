"""
Choosing a source for ∆T
========================

1. `ftp://ftp.iers.org/products/eop/rapid/standard/`_
2. `ftp://ftp.iers.org/products/eop/rapid/standard/csv/`_

* The CSV files are slightly larger, both compressed and uncompressed.
"""

#!/usr/bin/env python3

import argparse
import sys
from skyfield.api import load
from numpy import loadtxt

def main(argv):
    # parser = argparse.ArgumentParser(description=put description here)
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    # args = parser.parse_args(argv)
    # print(args.accumulate(args.integers))
    #for line in open('
    f = load.open('finals.all')
    import numpy as np
    import re
    r = re.compile(b'^(..)(..)(..)' + b'.' * 52 + b'(.\d........)', re.M)
    txt = f.read()
    print(r.findall(txt)[:10])
    print(r.findall(txt)[-10:])
    tups = r.findall(txt)
    for tup in tups:
        if any(field.isspace() for field in tup):
            print('============', tup)

    ts = load.timescale()

    # Old data.

    #t_old = ts.utc(1973, 1, range(365 * 30))

    # New data.

    f = load.open('finals.all')
    data = np.fromregex(f, r, [
        ('year', np.float64),
        ('month', np.float64),
        ('day', np.float64),
        ('ut1_minus_utc', np.float64),
    ])
    y = data['year']
    y += 1900
    y[y < 1973] += 100

    t = ts.utc(data['year'], data['month'], data['day'])
    y_new = data['ut1_minus_utc']
    y_old = t.dut1

    #loadtxt(f, usecols=(0, 1, 2))#, 10))
    # for line in f:

    # import numpy as np
    # t = np.arange(0.0, 2.0, 0.01)
    # s = 1 + np.sin(2 * np.pi * t)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    #ax.plot(t, s, label='label', linestyle='--')
    ax.plot(t.J, y_old, '-', label='deltat.dat')
    ax.plot(t.J, y_new, ',', label='finals.all')

    ax.grid()
    ax.set(xlabel='Year', title='UT1 minus UTC')
    #ax.set_aspect(aspect=1.0)
    #ax.axhline(1.0)
    #ax.axvline(1.0)
    plt.legend()
    fig.savefig('tmp.png')

    #     print(float(line[58:68]))
    #     break
"""
73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916  P    44.969     .500     2.839     .300   .143000   .137000   .8075000      .000      .000  \n\\
"""

if __name__ == '__main__':
    main(sys.argv[1:])

