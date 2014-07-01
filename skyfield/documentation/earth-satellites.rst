
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
(TLE) page of the Celestrak web site.

Once you have acquired the two-line orbital elements,
simply read them from a file or paste them directly into your script
to find out whether a satellite is above your local horizon:

.. testcode::

    text = """
    ISS (ZARYA)             
    1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082
    2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473
    """

    from skyfield.api import JulianDate, earth

    bluffton = earth.topos('40.8939 N', '83.8917 W')
    tup = (2014, 1, 21, 11, 18, 7)

    sat = earth.satellite(text)
    alt, az, distance = bluffton(utc=tup).observe(sat).altaz()

    print(alt)
    print(az)
    print(distance.km)

.. testoutput::

    13deg 49' 57.8"
    357deg 51' 21.8"
    1281.20477925

Observing a satellite directly returns an :ref:`apparent`
but without having to apply aberration and deflection,
because satellites travel with the Earth around the Sun
in the same relativistic frame of reference.

Propagation Errors
==================

If a satellite returns position and velocity vectors
with the values ``(nan, nan, nan)`` instead of real numbers,
it is because you are asking Skyfield about a date and time
for which the TLE elements cannot generate a position.

Earth satellite orbits are always in flux.
Even when a satellite is not firing its thrusters to adjust course,
it is constantly acted upon by effects
that include tidal forces, the Moon’s gravity,
and the (eventually fatal) drag of the Earth’s atmosphere.

The standard SGP4 theory of satellite motion that Skyfield uses
is a rough enough model of the near-Earth environment
that it can only predict a satellite position
to within an accuracy of a few kilometers.
This error grows larger as a TLE element set becomes several weeks old,
until its predictions are no longer meaningful.

After building a satellite object,
you can it ask for the date and time
of the *epoch* moment when its TLE element set’s predictions
are most accurate:

.. testcode::

    text = """
    GOCE                    
    1 34602U 09013A   13314.96046236  .14220718  20669-5  50412-4 0   930
    2 34602 096.5717 344.5256 0009826 296.2811 064.0942 16.58673376272979
    """
    sat = earth.satellite(text)
    print(sat.epoch.utc_jpl())

.. testoutput::

    A.D. 2013-Nov-09 23:03:03.9479 UT

Skyfield is willing to generate positions
for dates quite far from a satellite’s epoch,
even if they are not likely to be meaningful.
But it cannot generate a position
beyond the point where the elements stop making physical sense.
At that point, the satellite will start returning coordinates
that are all the special floating-point value ``nan`` “not-a-number.”

There is another, special case,
which is illustrated by the GOCE satellite we have just loaded up:
the orbital elements collapse
after the object enters the Earth’s atmosphere
and its orbit decays completely.
The elements given above are for the moment just a few hundred minutes
before the debris re-entered our atmosphere,
and so we can watch the SGP4 derivation collapse
once the satellite elements reach that age.

Starting 300 minutes after at the TLE element set’s epoch moment,
we can stride forward in 100 minute increments
by passing a the standard ``utc`` tuple
as described in :ref:`date-arrays`:

x.. testcode::

    year, month, day, hour, minute, second = sat.epoch.utc
    offsets = [100, 200, 300, 400, 500, 600]
    # offsets = range(-100, -90)
    jd = JulianDate(utc=(year, month, day, hour, minute-1000))
    geocentric = sat.gcrs(jd)

    x, y, z = geocentric.position.km
    print(jd.utc_jpl())
    print(x)

x.. testoutput::

    

x.. testcode::

    ra, dec, distance = geocentric.radec(epoch='date')
    print(ra.hms())

x.. testoutput::

    

x.. testcode::

    print(geocentric.sgp4_error)

x.. testoutput::

    

.. testsetup::

    import numpy as np
    np.set_printoptions(suppress=True, precision=4)

.. testcleanup::

    import numpy as np
    np.set_printoptions(precision=8)
