"""Parse data files from https://www.iers.org/ the IERS.

* ``finals.data`` is ``finals.all`` but starting with 1992 instead of 1973.

"""
import numpy as np
import re

FINALS_ALL_URL = 'ftp://ftp.iers.org/products/eop/rapid/standard/finals.all'
_DUT1_RE = re.compile(b'^......(.........) ' + b'.' * 42 + b'(.\d........)',
                      re.M)

def parse_dut1_from_finals_all(f):
    data = np.fromregex(f, _DUT1_RE, [
        ('mjd', np.float32),
        ('dut1', np.float32),
    ])
    return data['mjd'], data['dut1']
