"""Parse data files from the International Earth Rotation Service.

See:
https://datacenter.iers.org/eop.php
ftp://cddis.gsfc.nasa.gov/pub/products/iers/readme.finals2000A

"""
import numpy as np
from ..constants import DAY_S

inf = float('inf')

# This regular expression must remain a plain string; attempting to
# compile it triggers a bug in older NumPy versions like 1.14.3:
# https://github.com/skyfielders/python-skyfield/issues/372
_R = (b'(?m)^......(.........) . '
      b'(.\d.......)......... '
      b'(.\d.......).........  '
      b'.(.\d........)')

def parse_x_y_dut1_from_finals_all(f):
    return np.fromregex(f, _R, [
        ('mjd_utc', float),
        ('x_arcseconds', float),
        ('y_arcseconds', float),
        ('dut1', float),
    ])

def install_polar_motion_table(ts, finals_data):
    t = ts.utc(1858, 11, 17.0 + finals_data['mjd_utc'])
    ts.polar_motion_table = (
        t.tt,
        np.array(finals_data['x_arcseconds']),
        np.array(finals_data['y_arcseconds']),
    )

def _build_timescale_arrays(mjd_utc, dut1):
    big_jumps = np.diff(dut1) > 0.9
    leap_second_mask = np.concatenate([[False], big_jumps])
    tt_minus_utc = np.cumsum(leap_second_mask) + 32.184 + 12.0

    tt_jd = mjd_utc + tt_minus_utc / DAY_S + 2400000.5
    delta_t = tt_minus_utc - dut1
    delta_t_recent = np.array([tt_jd, delta_t])

    leap_dates = 2400000.5 + np.concatenate([
        [-inf], [41317.0, 41499.0, 41683.0], mjd_utc[leap_second_mask], [inf],
    ])
    leap_offsets = np.arange(8.0, len(leap_dates) + 8.0)
    leap_offsets[0] = leap_offsets[1] = 10.0

    return delta_t_recent, leap_dates, leap_offsets

# Compatibility with older Skyfield versions:

def parse_dut1_from_finals_all(f):
    data = parse_x_y_dut1_from_finals_all(f)
    return data['mjd_utc'], data['dut1']
