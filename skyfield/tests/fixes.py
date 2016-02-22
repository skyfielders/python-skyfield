"""Helpers for making Skyfield tests stable."""

from functools import partial

import datetime
import skyfield.api
import skyfield.timelib

def setup(utc=(2015, 10, 11, 10)):
    skyfield.api.load = partial(skyfield.api.load, verbose=False)
    skyfield.timelib.Timescale.utcnow = lambda self: datetime.datetime(*utc)

def teardown():
    skyfield.api.load = skyfield.iokit.load
    skyfield.timelib.Timescale.utcnow = datetime.datetime.utcnow
