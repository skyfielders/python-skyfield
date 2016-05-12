"""Various constants required by Skyfield."""

# Definitions.
AU_M = 149597870700             # per IAU 2012 Resolution B2
AU_KM = 149597870.700
ASEC360 = 1296000.0
DAY_S = 86400.0

# Angles.
ASEC2RAD = 4.848136811095359935899141e-6
DEG2RAD = 0.017453292519943296
RAD2DEG = 57.295779513082321
tau = 6.283185307179586476925287  # lower case, for symmetry with math.pi

# Physics.
C = 299792458.0

# Earth and its orbit.
ANGVEL = 7.2921150e-5
ERAD = 6378136.6
IERS_2010_INVERSE_EARTH_FLATTENING = 298.25642

# Heliocentric gravitational constant in meters^3 / second^2, from DE-405.
GS = 1.32712440017987e+20

# Time.
T0 = 2451545.0
B1950 = 2433282.4235

C_AUDAY = C * DAY_S / AU_M
