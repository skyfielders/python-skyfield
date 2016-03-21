"""Helpers for making Skyfield tests stable."""

import datetime
import skyfield.api
import skyfield.timelib

def setup(utc=(2015, 10, 11, 10)):
    skyfield.api.load = skyfield.iokit.Loader('.', verbose=False)
    skyfield.timelib.Timescale.utcnow = lambda self: datetime.datetime(*utc)

def teardown():
    skyfield.api.load = skyfield.iokit.Loader('.')
    skyfield.timelib.Timescale.utcnow = datetime.datetime.utcnow
