"""

"""
import os
from datetime import date
from skyfield import io
from skyfield import timescales

cache = io.Cache('.')
cache.npy_directory = os.path.dirname(__file__)
cache.npy_files = {

    timescales.download_leapseconds:
    (date(2014, 1, 1), ('leap_dates.npy', 'leap_offsets.npy')),

    }
