
==================
 Earth Satellites
==================

.. currentmodule:: skyfield

Skyfield is able to predict the positions of Earth satellites
from the Two-Line Element (TLE) files published
by organizations like `CelesTrak`_.
But there several limitations to be aware of
when using Skyfield to generate positions
for artificial satellites in Earth orbit:

.. _Celestrak: https://celestrak.com/

1. Do not expect perfect agreement between
   any two pieces of software that are trying to predict
   satellite positions from TLE data files.
   As Vallado, Crawford, and Hujsak observe in Appendix B of their
   crucial paper `Revisiting Spacetrack Report #3`_:

      “The maximum accuracy for a TLE is
      limited by the number of decimal places in each field.
      In general, TLE data is accurate to about a kilometer or so
      at epoch and it quickly degrades.”

.. _Revisiting Spacetrack Report #3:
    https://celestrak.com/publications/AIAA/2006-6753/

2. Satellite elements go rapidly out of date.
   You will want to pay attention to the “epoch” —
   the date on which an element set is most accurate —
   of every TLE element set you use.
   Elements are only useful for a week or two
   on either side of the epoch date,
   and for dates outside of that range
   you will want to download a fresh set of elements.

3. Given the low accuracy of TLE elements,
   there is no point in calling the usual Skyfield
   :meth:`~skyfield.positionlib.Barycentric.observe()` method
   that repeatedly re-computes an object’s position
   to account for the light-travel time to the observer.
   As we will see below,
   the difference is irrelevant for Earth satellites
   and not worth the added expense of re-computing the position
   several times in a row.

You can find satellite element sets at the
`NORAD Two-Line Element Sets <http://celestrak.com/NORAD/elements/>`_
page of the Celestrak web site.

Beware that the two-line element (TLE) format is very rigid.
The meaning of each character
is based on its exact offset from the beginning of the line.
You must download and use the element set’s text
without making any change to its whitespace.

Skyfield loader objects offer a :meth:`~skyfield.iokit.Loader.tle()`
method that can download and cache a file full of satellite elements
from a site like Celestrak.
A popular observing target for satellite observers
is the International Space Station,
which is listed in their ``stations.txt`` file:

.. testsetup::

    open('stations.txt', 'w').write("""\
    ISS (ZARYA)             
    1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082
    2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473
    """)

.. testcode::

    from skyfield.api import Topos, load

    stations_url = 'http://celestrak.com/NORAD/elements/stations.txt'
    sats = load.tle(stations_url)
    satellite = sats['ISS (ZARYA)']
    print(satellite)

.. testoutput::

    EarthSatellite 'ISS (ZARYA)' number=25544 epoch=2014-01-20T22:23:04Z

The value shown for the “epoch” is the all-important date
on which this set of elements is most accurate,
and before or after which they go rapidly out of date.
You can access this value as an attribute of the object
in case your program wants to check how old the elements are:

.. testcode::

    print(satellite.epoch)
    print(satellite.epoch.utc_jpl())

.. testoutput::

    <Time tt=2456678.433463>
    A.D. 2014-Jan-20 22:23:04.0004 UT

If the epoch is too far in the past,
you can provide :meth:`~skyfield.iokit.Loader.tle()`
with the ``reload`` option to force it to download new elements
even if the file is already on disk.

.. testcode::

   ts = load.timescale()
   now = ts.utc(2014, 1, 21, 11, 18, 7)

   days = now - satellite.epoch
   print('Difference: {:.3f} days'.format(days))
   if abs(days) > 14:
       sats = load.tle(stations_url, reload=True)
       satellite = sats['ISS (ZARYA)']

.. testoutput::

    Difference: 0.538 days



.. testcode::

    bluffton = Topos('40.8939 N', '83.8917 W')

    difference = satellite - bluffton
    position = difference.at(t)

    alt, az, distance = position.altaz()
    print(alt)
    print(az)
    print(distance.km)

    from skyfield.constants import AU_KM
    position.position.au[1] += 1/AU_KM
    alt, az, distance = position.altaz()
    print(alt)
    print(az)
    print(distance.km)

.. testoutput::

    13deg 50' 46.6"
    358deg 48' 55.9"
    1280.53654286

.. testcode::

    # OVERLY EXPENSIVE APPROACH - Compute both the satellite
    # and observer positions relative to the Solar System
    # barycenter ("ssb"), then call observe() to compensate
    # for light-travel time.

    de421 = load('de421.bsp')
    earth = de421['earth']
    ssb_bluffton = earth + bluffton
    ssb_satellite = earth + satellite
    position2 = ssb_bluffton.at(t).observe(ssb_satellite).apparent()

    # After all that work, how big is the difference, really?

    difference_km = (position2 - position).distance().km
    print('Difference between the two positions:')
    print('{0:.3f} km'.format(difference_km))

    difference_angle = position.separation_from(position2)
    print('Angle between the two positions in the sky:')
    print('{}'.format(difference_angle))

.. testoutput::

    Difference between the two positions:
    0.091 km
    Angle between the two positions in the sky:
    00deg 00' 05.0"

To find out whether the satellite is above your local horizon,
you will want to ask for its altitude and azimuth.
Negative altitudes lie below your horizon,
while positive altitude places the satellite above the horizon:


You can also ask for the position
to be expressed as right ascension and declination
relative to the fixed axes of the ICRS,
or else in dynamical coordinates
that are relative to the actual position
of the celestial equator and equinox on the date in question.
See :doc:`positions` to learn more about these possibilities:

.. testcode::

    ra, dec, distance = position.radec()  # ICRS/J2000

    print(ra)
    print(dec)

.. testoutput::

    01h 54m 36.45s
    +62deg 51' 50.6"

.. testcode::

    ra, dec, distance = position.radec(epoch='date')

    print(ra)
    print(dec)

.. testoutput::

    01h 55m 39.18s
    +62deg 55' 57.4"

The standard SGP4 theory of satellite motion that Skyfield uses
is a rough enough model of the near-Earth environment
that it can only predict a satellite position
to within an accuracy of a few kilometers.
This error grows larger as a TLE element set becomes several weeks old,
until its predictions are no longer meaningful.

Detecting Propagation Errors
============================

After building a satellite object,
you can examine the *epoch* date and time
when the TLE element set’s predictions are most accurate.
The ``epoch`` attribute is a :class:`Time`,
so it supports all of the standard Skyfield date methods:

.. testcode::

    from skyfield.api import EarthSatellite

    text = """
    GOCE                    
    1 34602U 09013A   13314.96046236  .14220718  20669-5  50412-4 0   930
    2 34602 096.5717 344.5256 0009826 296.2811 064.0942 16.58673376272979
    """
    lines = text.strip().splitlines()

    sat = EarthSatellite(lines[1], lines[2], lines[0])
    print(sat.epoch.utc_jpl())

.. testoutput::

    A.D. 2013-Nov-10 23:03:03.9479 UT

Skyfield is willing to generate positions
for dates quite far from a satellite’s epoch,
even if they are not likely to be meaningful.
But it cannot generate a position
beyond the point where the elements stop making physical sense.
At that point, the satellite will return a position and velocity
``(nan, nan, nan)`` where all of the quantities
are the special floating-point value ``nan`` which means *not-a-number*.

When a propagation error occurs and you get ``nan`` values,
you can examine the ``sgp4_error`` attribute of the returned position
to learn the error that the SGP4 propagator encountered.

We can take as an example the satellite elements above.
They are the last elements ever issued for GOCE,
about one day before the satellite re-entered the atmosphere
after an extended and successful mission.
Because of the steep decay of its orbit,
the elements are valid over an unusually short period —
from just before noon on Saturday to just past noon on Tuesday:

.. image:: _static/goce-decay.png

By asking for GOCE’s position just before or after this window,
we can learn about the propagation errors
that are limiting this TLE set’s predictions:

.. testcode::

    geocentric = sat.at(ts.utc(2013, 11, 9))
    print('Before:')
    print(geocentric.position.km)
    print(geocentric.message)

    geocentric = sat.at(ts.utc(2013, 11, 13))
    print('\nAfter:')
    print(geocentric.position.km)
    print(geocentric.message)

.. testoutput::

    Before:
    [ nan  nan  nan]
    mean eccentricity -0.001416 not within range 0.0 <= e < 1.0

    After:
    [ nan  nan  nan]
    mrt 0.997178 is less than 1.0 indicating the satellite has decayed

If you use a ``Time`` array to ask about an entire range of dates,
then ``sgp4_error`` will be a sequence filled in with ``None``
whenever the SGP4 propagator was successful
and otherwise recording the propagator error:

.. testcode::

    from pprint import pprint

    geocentric = sat.at(ts.utc(2013, 11, [9, 10, 11, 12, 13]))
    pprint(geocentric.message)

.. testoutput::

    ['mean eccentricity -0.001416 not within range 0.0 <= e < 1.0',
     None,
     None,
     None,
     'mrt 0.997178 is less than 1.0 indicating the satellite has decayed']
