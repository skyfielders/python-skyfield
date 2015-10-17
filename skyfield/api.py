"""Top-level objects and functions offered by the Skyfield library.

Importing this ``skyfield.api`` module causes Skyfield to load up the
default JPL planetary ephemeris ``de421`` and create planet objects like
``earth`` and ``mars`` that are ready for your use.

"""
import de421
from datetime import datetime
from math import pi
from .constants import tau
from .iokit import load
from .starlib import Star
from .timelib import JulianDate, T0, now, utc
from .toposlib import Topos
from .units import Angle
from .named_stars import NamedStar

__all__ = ['Angle', 'JulianDate', 'NamedStar', 'Star', 'Topos',
           'datetime', 'load', 'now', 'utc', 'T0', 'pi', 'tau']

def build_ephemeris():
    from .data.horizons import festoon_ephemeris
    from .jpllib import Ephemeris
    ephemeris = Ephemeris(de421)
    festoon_ephemeris(ephemeris)
    return ephemeris

ephemeris = build_ephemeris()
del build_ephemeris

earth = ephemeris.earth
mars = ephemeris.mars
jupiter = ephemeris.jupiter
