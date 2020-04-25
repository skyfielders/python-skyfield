"""Exceptions specific to the Skyfield library."""

from functools import wraps

class DeprecationError(Exception):
    """Explain that a Skyfield feature has been removed."""

class OutOfRangeTimeError(ValueError):
    """
    This exception is thrown is the given time is out of tange of the supported times.
    It has two extra attributes:

    - `first_valid_time`: the first supported time
    - `last_valid_time`: the last supported time
    - `out_of_range_times`: a list of booleans expressing which of the given times were out of range
    """

    def __init__(self, message, first_valid_time, last_valid_time, out_of_range_times):
        self.args = message,
        self.first_valid_time = first_valid_time
        self.last_valid_time = last_valid_time
        self.out_of_range_times = out_of_range_times

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
