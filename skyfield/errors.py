"""Exceptions specific to the Skyfield library."""

from functools import wraps

class DeprecationError(Exception):
    """Explain that a Skyfield feature has been removed."""

class EphemerisRangeError(ValueError):
    """An ephemeris has been asked about positions outside its time range.

    Attributes:

    - `start_time`, `end_time`: the range of times supported by the segment
    - `time_mask`: Boolean array where ``True`` marks out-of-range times
    - `segment`: the ephemeris segment that was asked for positions

    """
    def __init__(self, message, start_time, end_time, time_mask, segment):
        self.args = message,
        self.start_time = start_time
        self.end_time = end_time
        self.time_mask = time_mask
        self.segment = segment

def raise_error_for_deprecated_time_arguments(method):
    @wraps(method)
    def wrapper(self, jd=None, utc=None, tai=None, tt=None, tdb=None):
        if utc or tai or tt or tdb:
            raise DeprecationError("""method {0}() no longer takes timescale keyword arguments

If you need to quickly get an old Skyfield script working again, simply
downgrade to Skyfield version 0.6.1 using a command like:

        pip install skyfield==0.6.1

Skyfield used to let you skip building a Time object.  You were
allowed to provide a keyword argument like "utc", "tai", "tt", or "tdb"
to any method needing a date and time:

        position = earth.at(utc=(1980, 4, 20))    # the old way

But this forced Skyfield to maintain secret global copies of several
time scale data files, that need to be downloaded and kept up to date
for constructing Time objects to and from UTC.  Skyfield now makes
this collection of "Timescale" data files explicit.  You can create one
with "load.timescale()" and then build times using its methods, which
have the same names as the old keyword arguments:

        from skyfield.api import load
        ts = load.timescale()
        position = earth.at(ts.utc(1980, 4, 20))  # the new way\
""".format(method.__name__))

        return method(self, jd)
    return wrapper
