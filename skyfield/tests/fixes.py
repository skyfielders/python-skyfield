"""Helpers for making Skyfield tests stable."""

import datetime as dt
import skyfield.api
import skyfield.timelib

_real_datetime_class = dt.datetime

class datetime(dt.datetime):
    _fake_now = None
    @staticmethod
    def now():
        return datetime._fake_now

def setup(utc=(2015, 10, 11, 10)):
    # Importing Pandas raises a TypeError once our fake datetime is
    # installed.  So instead of waiting, let's import it now.
    __import__('pandas')

    fake_now = dt.datetime(*utc)
    skyfield.api.load = skyfield.iokit.Loader('.', verbose=False)
    skyfield.timelib.Timescale._utcnow = lambda self: fake_now
    datetime._fake_now = fake_now
    dt.datetime = datetime

def teardown():
    skyfield.api.load = skyfield.iokit.Loader('.')
    skyfield.timelib.Timescale._utcnow = dt.datetime.utcnow
    dt.datetime = _real_datetime_class
