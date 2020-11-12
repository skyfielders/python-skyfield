# import numpy as np
# from novas import compat as novas
# from skyfield import api, iokit
# from skyfield.constants import ASEC2RAD
# from skyfield.framelib import ICRS_to_J2000 as B
# from skyfield.functions import mxmxm, mxv, rot_x, rot_y, rot_z, tau

import numpy as np
from numpy import sin, cos
from novas import compat as novas
from skyfield import api
from skyfield.constants import ASEC2RAD, T0
from skyfield.functions import mxv, rot_x, rot_y, rot_z, tau

# xp_arcseconds = 0.0
# yp_arcseconds = 0.0
xp_arcseconds = 11.0
yp_arcseconds = 22.0

# class PolarTimescale(api.Timescale):
#     def M(self):
#         P = self.precession_matrix()
#         N = self.nutation_matrix()
#         return mxmxm(N, P, B)

# iokit.Timescale = PolarTimescale  # cheat: inject this alternate class

ts = api.load.timescale()
t = ts.utc(2020, 11, 11, 23, 37)

# radius = 6378.1366
# top = api.Topos(latitude_degrees=0.0, longitude_degrees=0.0, elevation_m=0.0,
#                 x=xp_arcseconds, y=yp_arcseconds)
# top_t = top.at(t)
# radius = top_t.distance().km

itrs_vector = [1.1, 1.2, 1.3]
v = itrs_vector

if yp_arcseconds or xp_arcseconds:
    sprime = -47.0e-6 * (t.tdb - T0) / 36525.0
    tiolon = -sprime * ASEC2RAD

    if 0:
        # Try it the NOVAS way.  More or less cut-and-pasted from novas.c:

        xpole = xp_arcseconds * ASEC2RAD
        ypole = yp_arcseconds * ASEC2RAD

        sinx = sin (xpole);
        cosx = cos (xpole);
        siny = sin (ypole);
        cosy = cos (ypole);
        sinl = sin (tiolon);
        cosl = cos (tiolon);

        xx =  cosx * cosl;
        yx =  sinx * siny * cosl + cosy * sinl;
        zx = -sinx * cosy * cosl + siny * sinl;
        xy = -cosx * sinl;
        yy = -sinx * siny * sinl + cosy * cosl;
        zy =  sinx * cosy * sinl + siny * cosl;
        xz =  sinx;
        yz = -cosx * siny;
        zz =  cosx * cosy;

        v = mxv(np.array([
            [xx, yx, zx],
            [xy, yy, zy],
            [xz, yz, zz],
        ]), v)

    else:
        # Try using NumPy matrix multiply.

        v = mxv(rot_x(-yp_arcseconds * ASEC2RAD), v)
        v = mxv(rot_y(-xp_arcseconds * ASEC2RAD), v)
        v = mxv(rot_z(-tiolon), v)

v = mxv(rot_z(t.gast / 24.0 * tau), v)
v = mxv(t.MT, v)
print(tuple(v))

EQUINOX = 1
FULL = 0
GCRS = 0
novas_km = novas.ter2cel(t.whole, t.ut1_fraction, t.delta_t,
                         xp_arcseconds, yp_arcseconds, itrs_vector,
                         method=EQUINOX, accuracy=FULL, option=GCRS)
print(novas_km)
print(v - novas_km)
