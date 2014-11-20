"""Top-level objects and functions offered by the Skyfield library.

Importing this ``skyfield.api`` module causes Skyfield to load up the
default JPL planetary ephemeris ``de421`` and create planet objects like
``earth`` and ``mars`` that are ready for your use.

"""
import de421
from datetime import datetime
from .starlib import Star
from .timelib import JulianDate, T0, now, utc
from .units import Angle

def build_ephemeris():
    from .data.horizons import festoon_ephemeris
    from .jpllib import Ephemeris
    ephemeris = Ephemeris(de421)
    festoon_ephemeris(ephemeris)
    return ephemeris

ephemeris = build_ephemeris()
del build_ephemeris

sun = ephemeris.sun
mercury = ephemeris.mercury
venus = ephemeris.venus
earth = ephemeris.earth
moon = ephemeris.moon
mars = ephemeris.mars
jupiter = ephemeris.jupiter
saturn = ephemeris.saturn
uranus = ephemeris.uranus
neptune = ephemeris.neptune
pluto = ephemeris.pluto

eight_planets = (mercury, venus, earth, mars, jupiter, saturn, uranus, neptune)
nine_planets = eight_planets + (pluto,)
