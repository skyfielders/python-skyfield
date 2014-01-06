import de421
from .jpllib import Ephemeris

ephemeris = Ephemeris(de421)
del Ephemeris

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

eight_planets = [mercury, venus, earth, mars, jupiter, saturn, uranus, neptune]
nine_planets = eight_planets + [pluto]
