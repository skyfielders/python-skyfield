
==============================
 Elegant Astronomy for Python
==============================

The Skyfield library provides a simple and safe interface
to high-precision astronomical routines
that are written in pure Python,
with NumPy as their only binary dependency;
this library itself requires no compilation.

You can install the current development version of Skyfield with::

    pip install https://github.com/brandon-rhodes/python-skyfield/archive/master.zip

The interface is still evolving,
but you can already ask the position of a planet in the sky
using one of three simple techniques::

    import skyfield.positionlib
    import skyfield.planets
    import de421

    ephemeris = skyfield.planets.Ephemeris(de421)
    earth = ephemeris.earth
    mars = ephemeris.mars

    # Getting the Julian date is not yet convenient.

    from sgp4.ext import jday
    t = jday(2013, 1, 30, 16, 8, 23)  # Julian date

    # What is the sky-catalog coordinate of Mars?

    ra, dec, d = mars.observe_from(earth(t)).radec()

    # Where is Mars in the sky, viewed from Earth's center?

    radec = mars.observe_from(earth(t)).apparent()

    # Where is Mars in the sky, viewed from Boston?

    boston = positionlib.Topos('71.0603 W', '42.3583 N',
        elevation=0.0, temperature=10.0, pressure=1010.0)
    radec = mars.observe_from(boston(t)).apparent()

All of these routines are designed
to operate very efficiently if,
instead of providing them with a single Julian date,
you provide them with a NumPy array of dates instead::

    import numpy as np

    t0 = jday(2013, 1, 30, 0, 0, 0)
    t = np.arange(t0, t0 + 365.0, 1.0)

    radec = mars.observe_from(boston(t)).apparent()

I will expand this README as the API continues to develop,
but I wanted to at least get these few details typed up
to help programmers who are already interested
in trying out Skyfield!
