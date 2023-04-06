#
# This is a quick investigation of the behavior of the geographic
# elevation formula `R / cos(lat) - aC`, which misbehaves at 90° north,
# as pointed out in GitHub issue #842.  As we 'add nines' to 89.0°, here
# is how the quantities behave:
#
# aC: this looks to be the 'local radius of curvature', of all things,
#     and it converges quickly to 6399593.625758493 m at 89.99999°, and
#     stays there, rock-solid, all the way to 90°.
#
# R / cos(lat): this is our problem, the value that is going sideways.
#     Converted to meters (* AU_M), it starts at a value very close to
#     aC, then begins to diverge - with `R / cos(lat) - aC` growing by
#     10× for each step from 89.99 through 89.999999999 before its
#     growth slows - but by that point the error in elevation is already
#     40 m.  So let's investigate the behavior of its components, then
#     return to this value below.
#
#  R: The radius of the line in the xy-plane connecting the pole with
#     the end of the distance vector.
#
# cos(lat): The scaling factor, at the angle `lat`, between the radius
#     of the unit circle, and a line drawn orthogonally from the pole's
#     line over to the tip of the position vector.
#
# R / cos(lat): So, we return to our ratio, now with the understanding
#     that it's trying to rebuild a kind-of-length for the position
#     vector, but using the computed `lat` instead of the simple true
#     angle between the position vector and the pole's vector.  So how
#     does this compare with the vector's true length?
#
# R / cos(lat) - length_of(xyz_au): This difference starts to converge
#     on 42841.31 before it starts to diverge at around 89.999999.  So
#     `R / cos(lat)` exceeds the actual vector length, because .
#
# Looking inside the _compute_latitude() routine, here are a few values:
#
# z:  The literal, unchanged z-coordinate of the position vector, giving
#     height above the xy-plane of the equator.
#
# aC * e2_sin_lat: The distance in the z-direction between the
#     position's actual z-coordinate, and the z-coordinate that the
#     position would have if it were on a sphere whose radius, the whole
#     way around, was the local ellipsoid's radius at the Earth latitude
#     under consideration.
#
# z + aC * e2_sin_lat: Per the previous two descriptions, this is the
#     position's z-coordinate adjusted 'upward' (north along the z-axis
#     for the northern hemisphere) to simulate the shape the vector
#     would have if the whole Earth had the radius of curvature of this
#     position on the ellipsoid.  And so it's this artificial z-value
#     that can be tossed into arctan2() with the literal xy-radius `R`.
#
# A few quick constants that will come up in various print() calls you
# might add to the code:
#
# 6399.6 km - Earth radius of curvature at pole
# 6378 km - Earth radius at equator
# 6357 km - Earth radius at pole
# 6335.4 km - Earth radius of curvature at equator

from numpy import arctan2, cos, sin, sqrt
from skyfield.api import load, wgs84, tau
from skyfield.constants import AU_M, DEG2RAD, RAD2DEG
from skyfield.functions import length_of

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
    xyz_au, x, y, aC, R, lat = self._compute_latitude(position)
    _, _, z = xyz_au

    # Official formula.
    # Jumps around, because aC is jumping around.
    height_au = R / cos(lat) - aC

    # So, an alternative?
    R2 = R
    lat2 = lat

    R2 = max(R, 11 / AU_M)
    lat2 = arctan2(z + aC * self._e2 * sin(lat2), R2)
    #lat2 = arctan2(z + aC * self._e2 * sin(lat2), R2)

    height_au = R2 / cos(lat2) - aC
    height_m = height_au * AU_M

    out = lat * RAD2DEG
    print(
        '{:<17s} {:<17s} {:20.12f}'.format(str(latitude), str(out), height_m),
        # aC * AU_M,  # Converges quickly on 6399 km, the radius of curvature.
        # R * AU_M,  # Dives towards zero.
        # cos(lat),  # Dives towards zero.
        (R / cos(lat)) * AU_M,  # Goes all wonky and jumps around.
        z * AU_M,  # Converges quickly on 6357 km, the actual polar radius.

        # Difference should be around 6399 - 6357 = ~42 km.
        # (R / cos(lat) - length_of(xyz_au)) * AU_M,
    )
