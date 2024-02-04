# -*- coding: utf-8 -*-

import numpy as np

import jplephem
import sgp4
import skyfield
from skyfield.api import load
from skyfield.functions import load_bundled_npy

def main():
    print('Skyfield version: {0}'.format(skyfield.__version__))
    print('jplephem version: {0}'.format(jplephem.__version__))
    print('sgp4 version: {0}'.format(sgp4.__version__))

    ts = load.timescale()
    fmt = '%Y-%m-%d'

    final_leap = (ts._leap_tai[-1] - 1) / (24 * 60 * 60)
    print('Built-in leap seconds table ends with leap second at: {0}'
          .format(ts.tai_jd(final_leap).utc_strftime()))

    arrays = load_bundled_npy('iers.npz')
    daily_tt = arrays['tt_jd_minus_arange']
    daily_tt += np.arange(len(daily_tt))
    start = ts.tt_jd(daily_tt[0])
    end = ts.tt_jd(daily_tt[-1])
    print('Built-in ∆T table from finals2000A.all covers: {0} to {1}'
          .format(start.utc_strftime(fmt), end.utc_strftime(fmt)))

main()
