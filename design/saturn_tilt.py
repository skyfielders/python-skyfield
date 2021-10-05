"""Compute the tilt of Saturn, for planetary_magnitude()."""

from skyfield.api import position_of_radec

# pc = PlanetaryConstants()
# pc.read_text(load('pck00008.tpc'))

# TODO: Figure out how to extract this from "pck00008.tpc".  For now,
# these values are copied and pasted.

BODY699_POLE_RA = 40.589
BODY699_POLE_DEC = 83.537

p = position_of_radec(40.589 / 15.0, 83.537, 1.0).xyz.au
print(p)
