# -*- coding: utf-8 -*-
"""Routines for interpreting data from the IAU Minor Planet Center."""

import io
import re

import pandas as pd

from ..data.spice import inertial_frames
from ..keplerlib import _KeplerOrbit
from ..timelib import julian_day

MPCORB_URL = 'https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz'

_MPCORB_COLUMNS = [
    ('designation_packed', (0, 7)),
    ('magnitude_H', (8, 13)),
    ('magnitude_G', (14, 19)),
    ('epoch_packed', (20, 25)),
    ('mean_anomaly_degrees', (26, 35)),
    ('argument_of_perihelion_degrees', (37, 46)),
    ('longitude_of_ascending_node_degrees', (48, 57)),
    ('inclination_degrees', (59, 68)),
    ('eccentricity', (70, 79)),
    ('mean_daily_motion_degrees', (80, 91)),
    ('semimajor_axis_au', (92, 103)),
    ('uncertainty', (105, 106)),
    ('reference', (107, 116)),
    ('observations', (117, 122)),
    ('oppositions', (123, 126)),
    ('observation_period', (127, 136)),
    ('rms_residual_arcseconds', (137, 141)),
    ('coarse_perturbers', (142, 145)),
    ('precise_perturbers', (146, 149)),
    ('computer_name', (150, 160)),
    ('hex_flags', (161, 165)),
    ('designation', (166, 194)),
    ('last_observation_date', (194, 202)),
]
_MPCORB_NECESSARY_COLUMNS = {
    'designation_packed', 'epoch_packed', 'mean_anomaly_degrees',
    'argument_of_perihelion_degrees', 'longitude_of_ascending_node_degrees',
    'inclination_degrees', 'eccentricity', 'mean_daily_motion_degrees',
    'semimajor_axis_au',
}
_MPCORB_DTYPES = {
    # These seem to be ignored by read_fwf()?
    'epoch_packed': 'category',
    'uncertainty': 'category',
    'coarse_perturbers': 'category',
    'precise_perturbers': 'category',
    'computer_name': 'category',
}
_MPCORB_CONVERTERS = {
    'designation_packed': str,
    'hex_flags': str,
}

def load_mpcorb_dataframe(fobj):
    """Parse a Minor Planet Center orbits file into a Pandas dataframe.

    See :doc:`kepler-orbits.rst`.  The MPCORB file format is documented
    at: https://minorplanetcenter.net/iau/info/MPOrbitFormat.html

    """
    # See https://github.com/pandas-dev/pandas/issues/18035
    fobj = io.TextIOWrapper(fobj)
    columns = _MPCORB_COLUMNS
    # if not slow:
    #     keepers = _MPCORB_NECESSARY_COLUMNS
    #     columns = [tup for tup in columns if tup[0] in keepers]
    names, colspecs = zip(*columns)
    df = pd.read_fwf(fobj, colspecs, names=names, dtypes=_MPCORB_DTYPES,
                     converters=_MPCORB_CONVERTERS)
    return df

def mpcorb_orbit(row, ts, gm_km3_s2):
    a = row.semimajor_axis_au
    e = row.eccentricity
    p = a * (1.0 - e*e)

    def n(c):
        return ord(c) - (48 if c.isdigit() else 55)

    def d(s):
        year = 100 * n(s[0]) + int(s[1:3])
        month = n(s[3])
        day = n(s[4])
        return julian_day(year, month, day) - 0.5

    epoch_jd = d(row.epoch_packed)
    t_epoch = ts.tt_jd(epoch_jd)

    minor_planet = _KeplerOrbit._from_mean_anomaly(
        p,
        e,
        row.inclination_degrees,
        row.longitude_of_ascending_node_degrees,
        row.argument_of_perihelion_degrees,
        row.mean_anomaly_degrees,
        t_epoch,
        gm_km3_s2,
        10,
        row.designation,
    )
    minor_planet._rotation = inertial_frames['ECLIPJ2000'].T
    return minor_planet

COMET_URL = 'https://www.minorplanetcenter.net/iau/MPCORB/CometEls.txt'

_COMET_COLUMNS = [
    ('number', (0, 4)),
    ('orbit_type', (4, 5)),
    ('designation_packed', (5, 12)),
    ('perihelion_year', (14, 18)),
    ('perihelion_month', (19, 21)),
    ('perihelion_day', (22, 29)),
    ('perihelion_distance_au', (30, 39)),
    ('eccentricity', (41, 49)),
    ('argument_of_perihelion_degrees', (51, 59)),
    ('longitude_of_ascending_node_degrees', (61, 69)),
    ('inclination_degrees', (71, 79)),
    ('perturbed_epoch_year', (81, 85)),
    ('perturbed_epoch_month', (85, 87)),
    ('perturbed_epoch_day', (87, 89)),
    ('magnitude_H', (91, 95)),
    ('magnitude_G', (96, 100)),
    ('designation', (102, 158)),
    ('reference', (159, 168)),
]
_COMET_FAST_COLUMNS = (
    'perihelion_year', 'perihelion_month', 'perihelion_day',
    'perihelion_distance_au', 'eccentricity', 'argument_of_perihelion_degrees',
    'longitude_of_ascending_node_degrees', 'inclination_degrees',
    'magnitude_H', 'magnitude_G',
    'designation',
)

_fast_comet_re = None
_fast_comet_sub = None

def load_comets_dataframe(fobj):
    """Parse a Minor Planet Center comets file into a Pandas dataframe.

    This imports only the fields essential for computing comet orbits.
    See :func:`~skyfield.data.mpc.load_comets_dataframe_slow()` for a
    slower routine that includes every comet data field.

    See :doc:`kepler-orbits.rst`.  The comet file format is documented
    at: https://www.minorplanetcenter.net/iau/info/CometOrbitFormat.html

    """
    global _fast_comet_re, _fast_comet_sub

    text = fobj.read()

    if _fast_comet_re is None:
        # Build a regular expression that will turn the fixed-width
        # comet file into a CSV that Pandas can import efficiently.

        keepers = set(_COMET_FAST_COLUMNS)
        pat = ['^']
        previous_end = None
        for name, (start, end) in _COMET_COLUMNS:
            if previous_end is not None:
                pat.append(' ' * (start - previous_end))
            keep = name in keepers
            if name == 'designation':
                pat.append('(.*?)  .*')
                break
            else:
                if keep:
                    pat.append('(')
                pat.append('.' * (end - start))
                if keep:
                    pat.append(')')
            previous_end = end

        pat = ''.join(pat)
        sub = ','.join('\\{}'.format(i + 1) for i in range(len(keepers)))

        _fast_comet_re = re.compile(pat.encode('ascii'), re.M)
        _fast_comet_sub = sub.encode('ascii')

    text = _fast_comet_re.sub(_fast_comet_sub, text)
    df = pd.read_csv(io.BytesIO(text), header=None, names=_COMET_FAST_COLUMNS)
    return df

def load_comets_dataframe_slow(fobj):
    """Parse a Minor Planet Center comets file into a Pandas dataframe.

    This routine reads in every single field from the comets data file.
    See :func:`~skyfield.data.mpc.load_comets_dataframe()` for a faster
    routine that omits some of the more expensive comet fields.

    See :doc:`kepler-orbits.rst`.  The comet file format is documented
    at: https://www.minorplanetcenter.net/iau/info/CometOrbitFormat.html

    """
    fobj = io.StringIO(fobj.read().decode('ascii'))
    names, colspecs = zip(*_COMET_COLUMNS)
    df = pd.read_fwf(fobj, colspecs, names=names)
    return df

def comet_orbit(row, ts, gm_km3_s2):
    e = row.eccentricity
    if e == 1.0:
        p = row.perihelion_distance_au * 2.0
    else:
        a = row.perihelion_distance_au / (1.0 - e)
        p = a * (1.0 - e*e)
    t_perihelion = ts.tt(row.perihelion_year, row.perihelion_month,
                         row.perihelion_day)

    comet = _KeplerOrbit._from_periapsis(
        p,
        e,
        row.inclination_degrees,
        row.longitude_of_ascending_node_degrees,
        row.argument_of_perihelion_degrees,
        t_perihelion,
        gm_km3_s2,
        10,
        row['designation'],
    )
    comet._rotation = inertial_frames['ECLIPJ2000'].T
    return comet

def unpack(designation_packed):
    def n(c):
        return ord(c) - (48 if c.isdigit() else 55)
    s = designation_packed
    s1 = s[1]
    if s1 == '/':
        return s
    return '{0[0]}/{1}{0[2]}{0[3]} {0[4]}{2}{3}'.format(
        s, n(s1), int(s[5:7]), s[7].lstrip('0'))
