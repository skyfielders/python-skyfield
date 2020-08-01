import locale
import numpy as np
from datetime import datetime, timedelta
from threading import Lock

from .timelib import Timescale, julian_date

_lock = Lock()

def _build_timescale(deltat_data, deltat_preds, leap_second_dat):
    data_end_time = deltat_data[0, -1]
    i = np.searchsorted(deltat_preds[0], data_end_time, side='right')
    delta_t_recent = np.concatenate([deltat_data, deltat_preds[:,i:]], axis=1)

    leap_dates, leap_offsets = leap_second_dat
    return Timescale(delta_t_recent, leap_dates, leap_offsets)

def parse_deltat_data(fileobj):
    """Parse the United States Naval Observatory ``deltat.data`` file.

    Each line file gives the date and the value of Delta T::

    2016  2  1  68.1577

    This function returns a 2xN array of raw Julian dates and matching
    Delta T values.

    """
    array = np.loadtxt(fileobj)
    year, month, day, delta_t = array.T
    return np.array((julian_date(year, month, day), delta_t))

def parse_deltat_preds(fileobj):
    """Parse the United States Naval Observatory ``deltat.preds`` file.

    The old format supplies a floating point year, the value of Delta T,
    and one or two other fields::

    2015.75      67.97               0.210         0.02

    The new format adds a modified Julian day as the first field:

    58484.000  2019.00   69.34      -0.152       0.117

    This function returns a 2xN array of raw Julian dates and matching
    Delta T values.

    """
    lines = iter(fileobj)
    header = next(lines)

    if header.startswith(b'YEAR'):
        # Format in use until 2019 February
        next(lines)             # discard blank line
        year_float, delta_t = np.loadtxt(lines, usecols=[0, 1]).T
    else:
        # Format in use since 2019 February
        year_float, delta_t = np.loadtxt(lines, usecols=[1, 2]).T

    year = year_float.astype(int)
    month = 1 + (year_float * 12.0).astype(int) % 12
    return np.array((julian_date(year, month, 1), delta_t))

def parse_leap_seconds(fileobj):
    """Parse the IERS file ``Leap_Second.dat``.

    The leap dates array can be searched with::

        index = np.searchsorted(leap_dates, jd, 'right')

    The resulting index allows (TAI - UTC) to be fetched with::

        offset = leap_offsets[index]

    """
    lines = iter(fileobj)
    for line in lines:
        if line.startswith(b'#  File expires on'):
            break
    else:
        raise ValueError('Leap_Second.dat is missing its expiration date')
    line = line.decode('ascii')

    with _lock:
        original_locale = locale.setlocale(locale.LC_ALL)
        locale.setlocale(locale.LC_ALL, 'C')
        try:
            dt = datetime.strptime(line, '#  File expires on %d %B %Y\n')
        finally:
            locale.setlocale(locale.LC_ALL, original_locale)

    # The file went out of date at the beginning of July 2016, and kept
    # downloading every time a user ran a Skyfield program.  So we now
    # build in a grace period:
    grace_period = timedelta(days=30)
    expiration_date = dt.date() + grace_period
    mjd, day, month, year, offsets = np.loadtxt(lines).T
    leap_dates = np.ndarray(len(mjd) + 2)
    leap_dates[0] = '-inf'
    leap_dates[1:-1] = mjd + 2400000.5
    leap_dates[-1] = 'inf'
    leap_offsets = np.ndarray(len(mjd) + 2)
    leap_offsets[0] = leap_offsets[1] = offsets[0]
    leap_offsets[2:] = offsets
    return expiration_date, (leap_dates, leap_offsets)
