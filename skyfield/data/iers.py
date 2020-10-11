"""Parse data files from https://www.iers.org/ the IERS.

* ``finals.data`` is ``finals.all`` but starting with 1992 instead of 1973.

"""
import numpy as np
import re

FINALS_ALL_URL = 'ftp://ftp.iers.org/products/eop/rapid/standard/finals.all'
_DUT1_RE = re.compile(b'^(..)(..)(..)' + b'.' * 52 + b'(.\d........)', re.M)

def parse_dut1_from_finals_all(f):
    data = np.fromregex(f, _DUT1_RE, [
        ('year', np.int16),
        ('month', np.int8),
        ('day', np.int8),
        ('dut1', np.float32),
    ])
    y = data['year']
    y += 1900
    y[y < 1973] += 100
    return y, data['month'], data['day'], data['dut1']
