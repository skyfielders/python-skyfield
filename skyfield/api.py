"""Top-level objects and functions offered by the Skyfield library.

Importing this library is not always the fastest way to use a Skyfield
feature, since importing this module involves importing almost the
entirety of Skyfield and its dependencies, but is the most convenient
way for most users to use Skyfield's main features.

"""
from datetime import datetime
from .constants import B1950, T0, pi, tau
from .constellationlib import load_constellation_map, load_constellation_names
from .iokit import Loader, load_file
from .planetarylib import PlanetaryConstants
from .positionlib import SSB, position_from_radec, position_of_radec
from .starlib import Star
from .sgp4lib import EarthSatellite
from .timelib import (
    GREGORIAN_START, GREGORIAN_START_ENGLAND, Time, Timescale, utc
)
from .toposlib import Topos, iers2010, wgs84
from .units import Angle, Distance, Velocity, wms

load = Loader('.')
N = E = +1.0
S = W = -1.0

__all__ = [
    'Angle', 'B1950', 'Distance', 'E', 'EarthSatellite',
    'GREGORIAN_START', 'GREGORIAN_START_ENGLAND',
    'Loader', 'PlanetaryConstants', 'N', 'S', 'SSB', 'Star', 'W',
    'T0', 'Time', 'Timescale', 'Topos', 'Velocity',
    'datetime', 'iers2010', 'load', 'load_constellation_map',
    'load_constellation_names', 'load_file',
    'position_from_radec', 'position_of_radec',
    'utc', 'pi', 'tau', 'wgs84', 'wms',
]
