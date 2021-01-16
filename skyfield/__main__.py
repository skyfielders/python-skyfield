# -*- coding: utf-8 -*-

import pkg_resources
import skyfield
from skyfield.api import load
from skyfield.functions import load_bundled_npy

def main():
    print('Skyfield version: {0}'.format(skyfield.__version__))
    print('jplephem version: {0}'.format(version_of('jplephem')))
    print('sgp4 version: {0}'.format(version_of('sgp4')))

    ts = load.timescale()
    fmt = '%Y-%m-%d'

    final_leap = (ts._leap_tai[-1] - 1) / (24 * 60 * 60)
    print('Built-in leap seconds table ends with leap second at: {0}'
          .format(ts.tai_jd(final_leap).utc_strftime()))

    arrays = load_bundled_npy('iers.npz')
    tt, delta_t = arrays['delta_t_recent']
    start = ts.tt_jd(tt[0])
    end = ts.tt_jd(tt[-1])
    print('Built-in ∆T table from finals2000A.all covers: {0} to {1}'
          .format(start.utc_strftime(fmt), end.utc_strftime(fmt)))

def version_of(distribution):
    try:
        d = pkg_resources.get_distribution(distribution)
    except pkg_resources.DistributionNotFound:
        return 'Unknown'
    else:
        return d.version

main()
