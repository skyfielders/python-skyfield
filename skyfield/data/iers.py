"""Parse data files from https://www.iers.org/"""

import numpy as np
import re

FINALS_ALL_URL = 'ftp://ftp.iers.org/products/eop/rapid/standard/finals.all'
_DUT1 = re.compile(b'^......(.........) ' + b'.' * 42 + b'(.\d........)', re.M)
inf = float('inf')

def parse_dut1_from_finals_all(f):
    data = np.fromregex(f, _DUT1, [
        ('mjd', np.float32),
        ('dut1', np.float32),
    ])
    return data['mjd'], data['dut1']

def build_timescale_arrays(mjd, dut1):
    big_jumps = np.diff(dut1) > 0.9
    leap_second_indices = np.concatenate([[False], big_jumps])
    delta_t = np.cumsum(leap_second_indices) - dut1 + 32.184 + 12.0
    delta_t_recent = np.array([mjd + 2400000.5, delta_t])

    leap_dates = 2400000.5 + np.concatenate([
        [-inf], [41317.0, 41499.0, 41683.0], mjd[leap_second_indices], [inf],
    ])
    leap_offsets = np.arange(8.0, len(leap_dates) + 8.0)
    leap_offsets[0] = leap_offsets[1] = 10.0

    return delta_t_recent, leap_dates, leap_offsets
