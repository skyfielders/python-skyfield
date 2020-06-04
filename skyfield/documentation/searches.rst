
==============================================
Searching for the dates of astronomical events
==============================================

As you plan your own observations or missions,
you might be interested in searching for the dates and times
of circumstances for which no ready-made solution
is already listed on Skyfield’s :doc:`almanac` page.

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

Finding discrete events
=======================

Perform a discrete event search
when you want to learn the date on which a continuous measurement,
like an angle or distance,
exceeds a particular value.

For example,
you might be reading a traditional astronomy text
and come across the definition of “quadrature” —
the moment when from the Earth’s point of view
a planet’s elongation from the Sun is 90°.

How might you compute the date of quadrature?

Always start by searching for the event yourself.
This will give you a sense for how the phenomenon behaves
before you attempt to unleash an automated search.
In the case of quadrature,
we might start by choosing Mars
and trying to compute its elongation today.
A dictionary, encyclopedia, or online reference
will clarify for us that a planet’s “elongation”
is its angular separation from the Sun,
so let’s compute the positions of Mars and the Sun
and then the angle between those positions:

.. testcode::

    from skyfield import api

    ts = api.load.timescale(builtin=True)
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
Maybe we came up
with a slightly different right ascension and declination for Mars.
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
which gives us high confidence that we are computing the elongation correctly.

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
It can help us develop an intuition
around how the value changes through time.

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
    plt.axes().axhline(90, color='r')  # Red line at 90°

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
and spends only a few months at a greater elongation.

Once we have learned to compute the value we are interested in
and have plotted its behavior,
there are three steps to solving for the dates on which it occurs:

1. Define a function of time returning an integer
   that changes each time the circumstance occurs.
   In a very simple case like this one,
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
    plt.tight_layout()
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
we are ready to unleash :func:`~skyfield.searchlib.find_discrete()`:

.. testcode::

    from skyfield.searchlib import find_discrete

    t1 = ts.utc(2018)
    t2 = ts.utc(2023)
    t, values = find_discrete(t1, t2, mars_quadrature)

    print(t)
    print(values)

.. testoutput::

    <Time tt=[2458202.1729387585 ... 2459818.728224116] len=5>
    [ True False  True False  True]

The result is a pair of arrays.
The first provides the dates and times of quadrature,
and the second provides the value
that our function switches to on each date.
The Python built-in function
`zip() <https://docs.python.org/3/library/functions.html#zip>`_
can iterate across both arrays at once
to pair up the dates with the values:

.. testcode::

    for ti, vi in zip(t, values):
        print(ti.utc_strftime('%Y-%m-%d %H:%M '), vi)

.. testoutput::

    2018-03-24 16:08  True
    2018-12-03 00:34  False
    2020-06-06 19:11  True
    2021-02-01 10:34  False
    2022-08-27 05:27  True

And we are done!
Those are the UTC dates
on which Mars reaches western quadrature
(when our discrete routine has just changed to ``True``)
and eastern quadrature
(when our routine has changed to ``False``),
as can be confirmed by comparing these dates
with those in a standard reference.

Finding extrema
===============

Sometimes you are not interested
in when a continuous function of time passes a threshold like 90°,
but when it reaches a minimum or maximum value —
the two possibilities are collectively called a function’s “extrema” —
whose exact value you might not be able to predict beforehand.

For example,
one challenge of observing Venus is that from Earth’s point of view
Venus’s smaller orbit
always keeps it within a few dozen degrees of the Sun.
Even when Venus is not so close to the Sun
that it’s hidden in the Sun’s glare,
it will be an evening star that’s already setting by the time we can see it
or a morning star that is soon followed by sunrise.

This leads observers to be interested in when Venus is farthest from the Sun —
when its elongation is greatest.

The steps are similar to those outlined in the previous section.
First, we define a function.

.. testcode::

    venus = eph['venus']

    def venus_elongation_degrees(t):
        e = earth.at(t)
        s = e.observe(sun).apparent()
        v = e.observe(venus).apparent()
        return s.separation_from(v).degrees

Then we compute a rough estimate
of how often Venus reaches greatest elongation.
The best approach is to generate a plot,
which will also give us a sense for how Venus’s elongation behaves.

.. testcode::

    from matplotlib import pyplot as plt

    plt.figure(figsize=(5, 2))
    plt.title('Elongation of Venus (degrees)')
    plt.xlabel('Year')
    plt.axes().grid(True)

    t = ts.utc(2018, 1, range(366 * 5))
    plt.plot(t.J, venus_elongation_degrees(t))

    plt.tight_layout()
    plt.savefig('venus-elongation.png')

.. image:: _static/venus-elongation.png

You might be surprised by the asymmetry between alternate minima —
between, say, the wide gradual minimum reached in mid-2019
versus the sharp quick minimum that comes next in mid-2020.
But if you investigate further,
in particular plotting Venus and the Earth in their orbits,
the reason will become clear:
Venus, on its faster orbit,
spends most of its time out on the other side of the Sun
gradually catching up with us,
then finally catches up and —
like a racecar zooming past us on the inside of a curve —
passes very quickly between our planet and the Sun,
generating the sharper “v” in our graph.

It looks like the maxima come no more often than each half-year,
so we can set the rough period to 180 days
and set the search routine to work.

.. testcode::

    from skyfield.searchlib import find_maxima

    venus_elongation_degrees.rough_period = 180.0

    t1 = ts.utc(2018)
    t2 = ts.utc(2023)
    t, values = find_maxima(t1, t2, venus_elongation_degrees)

    #print(t)
    #print(values)

..
   .. testoutput::

       <Time tt=[2458202.1729387585 ... 2459818.728224116] len=5>
       [ True False  True False  True]

.. testcleanup::

    import os
    os.rename('venus-elongation.png', '_static/venus-elongation.png')


.. testcode::

    for ti, vi in zip(t, values):
        print(ti.utc_strftime('%Y-%m-%d %H:%M:%S '), vi)

.. testoutput::

    2018-08-17 17:31:17  45.927744070690046
    2019-01-06 04:53:35  46.95601743723403
    2020-08-13 00:14:12  45.79100913313588
    2021-10-29 20:51:56  47.04509295323673
    2022-03-20 09:25:07  46.58628606386

(TODO)

Finding minima
--------------

(TODO)
