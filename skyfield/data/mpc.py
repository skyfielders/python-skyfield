# -*- coding: utf-8 -*-
"""Routines for interpreting data from the IAU Minor Planet Center."""

import pandas as pd

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

def load_mpcorb_dataframe(fobj, slow=False):
    """Parse a Minor Planet Center orbits file into a Pandas dataframe.

    The MPCORB file format is documented at:
    https://minorplanetcenter.net/iau/info/MPOrbitFormat.html

    """
    # TODO: should iokit handle decompression in open()?
    columns = _MPCORB_COLUMNS
    if not slow:
        columns = [tup for tup in columns
                   if tup[0] in _MPCORB_NECESSARY_COLUMNS]
    names, colspecs = zip(*columns)
    df = pd.read_fwf(fobj, colspecs, names=names, dtypes=_MPCORB_DTYPES,
                     converters=_MPCORB_CONVERTERS)
    #skiprows=43)
    return df

COMET_URL = 'https://www.minorplanetcenter.net/iau/MPCORB/CometEls.txt'

_COMET_SLOW_COLUMNS = [
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
_COMET_FAST_COLUMN_NAMES, _COMET_FAST_COLUMN_NUMBERS = zip(
    ('designation_packed', 0),
    ('perihelion_year', 1),
    ('perihelion_month', 2),
    ('perihelion_day', 3),
    ('perihelion_distance_au', 4),
    ('eccentricity', 5),
    ('argument_of_perihelion_degrees', 6),
    ('longitude_of_ascending_node_degrees', 7),
    ('inclination_degrees', 8),
    ('magnitude_H', 10),
    ('magnitude_G', 11),
)
_COMET_DTYPES = {
    'number': 'float',  # since older Pandas does not support NaN for integers
    'orbit_type': 'category',
    'perihelion_year': 'int',
}

import re
_pat = re.compile(br'^([^ ]*) +([a-z])', flags=re.M)
import io

def load_comets_dataframe(fobj):
    """Parse a Minor Planet Center comets file into a Pandas dataframe.

    The comet file format is documented at:
    https://www.minorplanetcenter.net/iau/info/CometOrbitFormat.html

    This uses a fast Pandas import routine on only the data fields
    essential for computing comet orbits, speeding up the import by a
    factor of 2 or 3.  But in return, each comet’s full name will be
    missing; only its packed designation is included as an identifier.

    See :func:`~skyfield.data.mpc.load_comets_dataframe_slow()` for a
    slower routine that includes every comet data field.

    """
    text = fobj.read()
    text = _pat.sub(br'\1\2', text)
    df = pd.read_csv(
        io.BytesIO(text), sep=r'\s+', header=None,
        names=_COMET_FAST_COLUMN_NAMES,
        usecols=_COMET_FAST_COLUMN_NUMBERS,
    )
    return df

def load_comets_dataframe_slow(fobj):
    """Parse a Minor Planet Center comets file into a Pandas dataframe.

    The comet file format is documented at:
    https://www.minorplanetcenter.net/iau/info/CometOrbitFormat.html

    This routine reads in every single field from the comets data file.
    See :func:`~skyfield.data.mpc.load_comets_dataframe()` for a faster
    routine that omits some of the more expensive comet fields.

    """
    fobj = io.StringIO(fobj.read().decode('ascii'))
    names, colspecs = zip(*_COMET_SLOW_COLUMNS)
    df = pd.read_fwf(fobj, colspecs, names=names, dtypes=_COMET_DTYPES)
    return df
