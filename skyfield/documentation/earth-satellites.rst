
==================
 Earth Satellites
==================

.. currentmodule:: skyfield

The orbital elements for human-launched Earth satellites
can go out of date in only a few days,
so it is best to start a programming or analysis session
by downloading fresh orbital elements for your object of interest
from the
`NORAD Two-Line Element Sets <http://celestrak.com/NORAD/elements/>`_
page of the Celestrak web site.

Once you have acquired the two lines of orbital elements,
simply read them from a file or paste them directly into your script
to find out whether a satellite is above your local horizon:

.. testcode::

    text = """
    ISS (ZARYA)             
    1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082
    2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473
    """

    from skyfield.api import earth

    bluffton = earth.topos('40.8939 N', '83.8917 W')
    tup = (2014, 1, 21, 11, 18, 7)

    sat = earth.satellite(text)
    alt, az, distance = bluffton(utc=tup).observe(sat).altaz()

    print(alt)
    print(az)
    print(distance.km)

.. testoutput::

    13deg 31' 26.0"
    21deg 34' 12.0"
    1296.58733908

Observing a satellite directly returns an :ref:`apparent`
but without having to apply aberration and deflection,
because satellites travel with the Earth around the Sun
in the same relativistic frame of reference.
