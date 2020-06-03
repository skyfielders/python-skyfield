
Searching for your own events
=============================

As you plan your own observations or missions,
you might be interested in searching for the dates and times
of circumstances for which there is no existing almanac routine.

There are two kinds of searches you can do.
One is an “event” search
that looks for a moment of position or alignment
that fits some definition that you are looking for.
The other is an “extremum” search
where you want to know when a value reaches its maximum or minimum value.
For the moment,
this documentation page only discusses searching for events;
check out the ``find_maxima()`` routine in the ``searchlib.py`` source code
for how Skyfield searches for maximum or minimum values.

Finding events
==============

You may sometimes be interested in astronomical circumstances
for which Skyfield’s :doc:`almanac` page
provides no ready-made solution.

For example,
you might be reading a traditional astronomy text
and come across the definition of “quadrature” —
the moment when from the Earth’s point of view
a planet’s elongation from the Sun is 90°.
How might you compute the date of quadrature yourself?

Always start by searching for the event yourself.
This will give you a sense for how the phenomenon behaves
before you attempt to unleash an automated search.
In the case of quadrature,
we might start by choosing Mars
and trying to compute its elongation today.
A dictionary, encyclopedia, or online reference
will clarify for us that a planet’s “elongation”
is its angular separation from the Sun,
so let’s compute their positions and then the angle between them:

.. testcode::

    from skyfield import api

    ts = api.load.timescale()
    t = ts.utc(2020, 6, 2)

    eph = api.load('de421.bsp')
    earth, sun, mars = eph['earth'], eph['sun'], eph['mars']

    e = earth.at(t)
    s = e.observe(sun)
    m = e.observe(mars)

    print('%.4f' % s.separation_from(m).degrees)

.. testoutput::

    88.5752

Have we computed the elongation correctly?
We should always double-check our work against other authorities when possible.
Given that same date,
the `NASA JPL HORIZONS <https://ssd.jpl.nasa.gov/horizons.cgi>`_ site
can produce a table with a “S-O-T” column,
which (as it explains in the definitions below the table)
is the solar elongation::

  Date__(UT)__HR:MN   R.A._____(ICRF)_____DEC    S-O-T /r
 ********************************************************
  2020-Jun-02 00:00   23 01 20.93 -08 51 51.6  88.5698 /L

Well, drat.

The HORIZONS system gives “88.5698” as the elongation at midnight UTC.
That’s not the same number we computed,
though it’s close.
How could Skyfield and HORIZONS have come up with different numbers
for the elongation?
Maybe they came up with a slightly different right ascension and declination.
Let’s check:

.. testcode::

    ra, dec, distance = m.radec()
    print(ra, '/', dec)

.. testoutput::

    23h 01m 20.93s / -08deg 51' 51.6"

Nope, that’s not the difference —
these Skyfield numbers are an exact match
for the right ascension and declination in the HORIZONS output shown above.

It is always worthwhile to study every detail of the HORIZONS output
when investigating a difference in result.
Here, for example, is the first paragraph of its definition
of that ``S-O-T`` field::

  S-O-T /r =
     Sun-Observer-Target angle; target's apparent SOLAR ELONGATION seen from
 the observer location at print-time. Angular units: DEGREES

So that’s the difference!
We computed the angle between the *astrometric* positions of the Sun and Mars,
whereas the elongation is more properly an angular difference
between *apparent* positions.
Thus:

.. testcode::

    s = e.observe(sun).apparent()
    m = e.observe(mars).apparent()

    print('%.4f' % s.separation_from(m).degrees)

.. testoutput::

    88.5698

Much better!
We now have a perfect match with HORIZONS
which gives us confidence that we are computing elongation correctly.

Next let’s search for a moment of quadrature.
I did not deliberately plan the example this way,
but it looks like Mars is very nearly near quadrature as I type this!
To determine whether quadrature was just reached
or is a few days in the future,
let’s compute the value over a few days
and see whether it’s growing or shrinking:

.. testcode::

    def mars_elongation_degrees(t):
        e = earth.at(t)
        s = e.observe(sun).apparent()
        m = e.observe(mars).apparent()
        return s.separation_from(m).degrees

    t = ts.utc(2020, 6, range(2 - 3, 2 + 3))

    for ti, ei in zip(t, mars_elongation_degrees(t)):
        print('%s %.4f' % (ti.utc_strftime('%b %d'), ei))

.. testoutput::

    May 30 87.6881
    May 31 87.9810
    Jun 01 88.2749
    Jun 02 88.5698
    Jun 03 88.8657
    Jun 04 89.1626

We see that the elongation of Mars is growing slowly right now,
at a rate of less than a degree per day,
but is very nearly at our target value of 90°.
Does is always grow slowly?
Does it wane at the same rate?
Are there periods during which its change is quick
and others during which it is slow?

I always recommend plotting any value
on which you are planning to perform a search.
It helps develop your intuition
around how the value changes.

.. testsetup::

    import matplotlib
    matplotlib.use('Agg')  # to avoid “no display name” error on Travis CI
    del matplotlib

.. testcode::

    from matplotlib import pyplot as plt

    plt.figure(figsize=(5, 3))
    plt.title('Elongation of Mars (degrees)')
    plt.xlabel('Year')
    plt.axes().grid(True)
    plt.axes().axhline(90, color='r')

    t = ts.utc(2018, 1, range(366 * 5))
    plt.plot(t.J, mars_elongation_degrees(t))

    plt.tight_layout()
    plt.savefig('mars-elongation.png')

.. image:: _static/mars-elongation.png

.. testcleanup::

    import os
    os.rename('mars-elongation.png', '_static/mars-elongation.png')

The dates of quadrature are where the elongation
intersects the red 90° line that we have drawn across the figure.
Mars seems to spend most of its time
with an elongation of less than 90° —
over on the same side of the sky as the Sun —
and then spend only a few months at a greater elongation.

Once we have learned to compute the value we are interested in
and have plotted its behavior,
there are three steps to solving for the dates on which it occurs:

1. Define a function of time returning an integer
   that changes each time the circumstance occurs.
   In a very simple case,
   you can simply use the values ``False`` and ``True``
   because in Python those are the integers zero and one.

2. Give the function a ``rough_period`` attribute
   telling the search routine how far apart to space its test dates
   when it first searches for where your function switches values.

3. Pass the function
   to the same :func:`~skyfield.searchlib.find_discrete()` routine
   that you would use for a search with the standard almanac functions.

The first step is quite easy in this case.
We simply need to compare the elongation with 90°.
This transforms the continuous angle measurement
into a discrete function
that jumps instantly between zero and one.

.. testcode::

    def mars_quadrature(t):
        e = earth.at(t)
        s = e.observe(sun).apparent()
        m = e.observe(mars).apparent()
        return s.separation_from(m).degrees >= 90

Since the Python values ``False`` and ``True``
are really the integers 0 and 1,
a plot of this function shows a square wave
whose positive excursions
identify the periods of time during which Mars is more than 90° from the Sun —
as we can verify by comparing this plot with our earlier plot.

.. testcode::

    from matplotlib import pyplot as plt

    plt.figure(figsize=(5, 1.5))
    plt.plot(t.J, mars_quadrature(t))
    plt.savefig('mars-quadrature.png')

.. image:: _static/mars-quadrature.png

.. testcleanup::

    import os
    os.rename('mars-quadrature.png', '_static/mars-quadrature.png')

The second step is to identify the ``rough_period``
over which our phenomenon cycles between true and false,
as measured in days.
Our plots suggest that the Mars elongation cycle
takes more than 2 years but less than 3 years.
Let’s use a round guess of 700 days.

.. testcode::

    mars_quadrature.rough_period = 700.0

Finally,
we are ready to unleash :func:`~skyfield.searchlib.find_discrete()`!

.. testcode::

    from skyfield.searchlib import find_discrete

    t1 = ts.utc(2018)
    t2 = ts.utc(2023)
    t, values = find_discrete(t1, t2, mars_quadrature)

    for ti, vi in zip(t, values):
        print(ti.utc_strftime('%Y-%m-%d '), vi)

.. testoutput::

    2018-03-24  True
    2018-12-03  False
    2020-06-06  True
    2021-02-01  False
    2022-08-27  True

And we are done!
Those are the UTC dates
on which Mars reaches western quadrature
(when our discrete routine has just changed to ``True``)
and eastern quadrature
(when our routine has changed to ``False``),
as can be confirmed by comparing these dates
with those in a standard reference.
