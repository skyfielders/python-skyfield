# -*- coding: utf-8 -*-
"""Routines to download Earth orientation data."""

import numpy as np
import os
from skyfield.api import load
from skyfield.timelib import julian_date

def morrison_and_stephenson_2004_table():
    """Table of smoothed Delta T values from Morrison and Stephenson, 2004."""
    import pandas as pd
    f = load.open('http://eclipse.gsfc.nasa.gov/SEcat5/deltat.html')
    tables = pd.read_html(f.read())
    df = tables[0]
    return pd.DataFrame({'year': df[0], 'delta_t': df[1]})

def parse_S15_table(f):
    """Parse polynomial coefficients from Table S15.

    The table is available at the website of Her Majesty's Nautical
    Almanac Office, from the paper “Measurement of the Earth's Rotation:
    720 BC to AD 2015” by L.V. Morrison, F.R. Stephenson, C.Y. Hohenkerk
    and M. Zawilski 2021.

    """
    # http://astro.ukho.gov.uk/nao/lvm/Table-S15.2020.txt
    content = f.read()
    banner = b'- ' * 36 + b'-\n'
    sections = content.split(banner)
    names = sections[1].splitlines()[-1].decode('utf-8').split()
    table = np.loadtxt(sections[2].splitlines())
    return names, table.T

def main():
    thisdir = os.path.dirname(__file__)

    df = morrison_and_stephenson_2004_table()
    year = df.year.values
    jd = julian_date(year, 1, 1)
    delta_t = df.delta_t.values
    array = np.array((jd, delta_t))
    np.save(os.path.join(thisdir, 'morrison_stephenson_deltat.npy'),
            array)

if __name__ == '__main__':
    main()
