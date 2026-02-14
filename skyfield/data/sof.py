# -*- coding: utf-8 -*-
"""Routines for interpreting data from Project Pluto standard orbit file (.sof)."""

import io
import re
import math

import pandas as pd

from ..data.spice import inertial_frames
from ..keplerlib import _KeplerOrbit

# example file:
# Name        |Tp      .       |Te      |q          |i  .      |Om .      |om .      |e         |rms |n_o  |Tfirst  |Tlast   |Perts  |H .  |G . ^
#        65803 20221021.7988553 20221030  1.01292656   3.413508  72.998191 319.555213 0.38357671 0.73  5582 19950307 20250227 M-v 3Ek 18.11  0.15

_SOF_COLUMNS = [
    ('name', (0, 12)),
    ('perihelion_date', (13, 29)),
    ('perturbed_epoch_date', (30, 38)),
    ('perihelion_distance_au', (39, 50)),
    ('inclination_degrees', (51, 61)),
    ('longitude_of_ascending_node_degrees', (62, 72)),
    ('argument_of_perihelion_degrees', (73, 84)),
    ('eccentricity', (84, 94)),
    ('rms', (95, 99)),
    ('number_of_orbits', (100, 105)),
    ('first_used_obs_date', (105, 114)),
    ('last_used_obs_date', (115, 123)),
    ('perturbators', (124, 131)),
    ('absolute_magnitude', (132, 137)),
    ('magnitude_slope', (138, 143))
]

_SOF_DTYPES = {
    'name': str, 'perihelion_date': str, 'perturbed_epoch_date': str,
    'first_used_obs_date': str, 'last_used_obs_date': str, 'perturbators': str
}

def load_sof_dataframe(fobj):
    """Parse a Project Pluto standard orbit file into a Pandas dataframe.

    This routine reads in every single field from the data file.
    See :func:`~skyfield.data.sof.load_sof_dataframe()` for a faster
    routine that omits some of the more expensive fields.

    See :doc:`kepler-orbits`.  The sof file format is documented at:
    https://projectpluto.com/orb_form.htm
    """
    fobj = io.StringIO(fobj.read().decode('ascii'))
    names, colspecs = zip(*_SOF_COLUMNS)
    df = pd.read_fwf(fobj, colspecs=colspecs, names=names, dtype=_SOF_DTYPES)
    return df

def sof_orbit(row, ts, gm_km3_s2):
    e = row.eccentricity
    if e == 1.0:
        p = row.perihelion_distance_au * 2.0
    else:
        a = row.perihelion_distance_au / (1.0 - e)
        p = a * (1.0 - e*e)
    # example date format: 20221021.7988553
    year = int(row.perihelion_date[0:4])
    month = int(row.perihelion_date[4:6])
    day = int(row.perihelion_date[6:8])
    day_fraction = float(row.perihelion_date[8:16])
    hour_fraction = day_fraction * 24
    hour = math.floor(hour_fraction)
    minute_fraction = (hour_fraction - hour) * 60
    minute = math.floor(minute_fraction)
    second = (minute_fraction - minute) * 60
    t_perihelion = ts.tt(year, month, day, hour, minute, second)

    sof = _KeplerOrbit._from_periapsis(
        p,
        e,
        row.inclination_degrees,
        row.longitude_of_ascending_node_degrees,
        row.argument_of_perihelion_degrees,
        t_perihelion,
        gm_km3_s2,
        10,
        row['name'],
    )
    sof._rotation = inertial_frames['ECLIPJ2000'].T
    return sof
