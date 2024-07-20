#
# This is a quick investigation of the behavior of the geographic
# elevation formula `R / cos(lat) - aC`, which misbehaves at 90° north,
# as pointed out in GitHub issue #842.  As we 'add nines' to 89.0°, here
# is how the quantities behave:
#
# aC: this is the 'local radius of curvature', which converges quickly
#     to 6399593.625758493 m at 89.99999°, and stays there, rock-solid,
#     all the way to 90°.
#
#  R: The radius of the line in the xy-plane connecting the pole with
#     the end of the distance vector, which (obviously) becomes very
#     very small as we approach the pole, though machine precision seems
#     to prevent its ever becoming zero.
#
# cos(lat): Alas, fairly imprecise, as the only way it can tell exactly
#     how far we are from the pole is to look at the very last digits of
#     angles like 89.9999999, which only have a few bits left beneath
#     them to hold their precision.
#
# R / cos(lat): So this is our problem, the value that is going
#     sideways.  For the Earth's surface at the North Pole, converted to
#     meters (* AU_M), it starts at a value very close to aC, then
#     begins to diverge - with `R / cos(lat) - aC` growing by 10× for
#     each step from 89.99 through 89.999999999 before its growth slows
#     - but by that point the error in elevation is already 40 m.
#
# How can we compute the height without recourse to the cosine of the
# latitude?  Drawing a diagram, we have a triangle with one leg that's
# the axis through the poles; one leg that's R, sticking out from the
# axis into the xy-plane; and a hypotenuse that's the radius of
# curvature plus the height above-ground-level of our target.
#
# R / cos(lat) is one way of computing the hypotenuse, yes.
#
# But we can also just use Pythagoras.  We already know one leg is R.
# The other?  It turns out that it's the quantity passed to arctan2() to
# compute the latitude!  Thus, `z + aC * e2_sin_lat`, which we can give
# a name to and pass as an additional return value.
#
# Some common values:
#
# 6399.6 km - Earth radius of curvature at pole
# 6378 km - Earth radius at equator
# 6357 km - Earth radius at pole
# 6335.4 km - Earth radius of curvature at equator

from numpy import cos, sqrt
from skyfield.api import load, wgs84
from skyfield.constants import AU_M, RAD2DEG

ts = load.timescale()
t = ts.utc(2023, 3, 2, 4, 13, 0)
self = wgs84  # 'self', as though inside method

for nines in range(16):
    #
    # Code below is adapted from geographic_position_of().
    #
    latitude = 90.0 - 1.0 / 10.0 ** nines
    # latitude = 1.0 - 1.0 / 10.0 ** nines  # To test behavior at equator.

    position = wgs84.latlon(latitude, 0.0, 0.0).at(t)
    xyz_au, x, y, R, aC, hyp, lat = self._compute_latitude(position)
    _, _, z = xyz_au

    # Official formula: jumps around, because R / cos(lat) is jumping
    # around, as both R and cos(lat) are becoming very small.
    height_au = R / cos(lat) - aC
    height_m = height_au * AU_M

    # So, new formula, building the hypotenuse with Pythagoras.
    thing = sqrt(hyp * hyp + R * R)
    height2_au = thing - aC
    height2_m = height2_au * AU_M

    out = lat * RAD2DEG
    print('{:<17s} {:<17s} {:20.12f} {:20.12f}'.format(
        str(latitude),
        str(out),
        height_m,
        height2_m,
    ))
