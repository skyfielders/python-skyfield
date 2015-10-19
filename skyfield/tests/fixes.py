"""Helpers for making Skyfield tests stable."""

from functools import partial

import skyfield.api

def setup(utc=(2015, 10, 11, 10)):
    def now():
        """Return a constant "now"."""
        return skyfield.api.JulianDate(utc=utc)
    skyfield.api.now = now
    skyfield.api.load = partial(skyfield.api.load, verbose=False)

def teardown():
    skyfield.api.now = skyfield.timelib.now
    skyfield.api.load = skyfield.iokit.load
