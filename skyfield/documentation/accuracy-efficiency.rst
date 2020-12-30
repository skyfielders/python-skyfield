
====================================
 Skyfield’s Accuracy and Efficiency
====================================

This document is a work in progress,
that will be expanded into a full guide.
Right now it covers only one topic.

.. _polar motion:

------------
Polar Motion
------------

It was discovered more than a century ago
that the Earth’s crust has a slight freedom
to wobble with respect to the Earth’s axis of rotation in space,
because the continents and ocean basis are bound to the planet’s mass
only through the fluid coupling of our planet’s viscous mantle.
In Skyfield you can see the size of the effect
by loading an official data file
from the International Earth Rotation Service (IERS)
and measuring the maximum excursions
of the polar motion parameters 𝑥 and 𝑦:

.. testcode::

    from skyfield.api import load
    from skyfield.data import iers

    url = load.build_url('finals2000A.all')
    with load.open(url) as f:
        finals_data = iers.parse_x_y_dut1_from_finals_all(f)

    ts = load.timescale()
    iers.install_polar_motion_table(ts, finals_data)

    tt, x_arcseconds, y_arcseconds = ts.polar_motion_table
    print('x:', max(abs(x_arcseconds)), 'arcseconds')
    print('y:', max(abs(y_arcseconds)), 'arcseconds')

.. testoutput::

    x: 0.32548 arcseconds
    y: 0.596732 arcseconds

In what kinds of Skyfield calculations
does the exact position of the Earth’s crust come into play?

* Polar motion affects the position of an observer on the Earth’s surface.

* Polar motion therefore also affects the relative position
  of a target with respect to an observer on the Earth’s surface.

* Polar motion directly affects altazimuth coordinates,
  since the polar angles 𝑥 and 𝑦 tilt the zenith and local horizon
  against which altitude and azimuth are measured for a particular observer.

* Finally,
  polar motion affects the position of any observation target
  that’s located on the Earth’s surface —
  for example, if you are calculating the position of a ground station
  from the perspective of a space probe.

To have Skyfield apply polar motion when computing positions and coordinates,
simply install the IERS tables on your timescale object
as shown in the example code above.
Polar motion will be used everywhere that it applies.
