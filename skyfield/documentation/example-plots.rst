
===============
 Example Plots
===============

This section of the documentation
will gradually accumulate example scripts
for producing images from Skyfield computations.

Note that these example scripts
are written for fairly recent versions of `matplotlib`_.
If you try running them on a system
with an older version of the library,
you might see errors — in particular with how they specify colors,
in which case you can try omitting those parameters
to get the script running.
In any case,
these are only intended to be a starting point
for building your own scripts,
either with matplotlib or whatever other plotting library you prefer.

.. _neowise-chart:

Drawing a finder chart for comet NEOWISE
========================================

Here is a stand-alone script
that brings together four different data sources —
a planetary ephemeris, a comet orbit database, a large star catalog,
and constellation diagrams —
to plot the course of Comet NEOWISE across Ursa Major
over one week of July 2020:

.. image:: _static/neowise-finder-chart.png

.. testsetup::

    import matplotlib
    matplotlib.use('Agg')  # to avoid “no display name” error on Travis CI
    del matplotlib

    import sys
    sys.path[0:0] = ['../../examples']
    import comet_neowise_chart

Its code includes many design decisions and presentation tweaks
that you will probably want to adjust for your own project.
Use the script as a starting point:

.. include:: ../../examples/comet_neowise_chart.py
   :literal:

.. testcleanup::

    import os
    os.rename('neowise-finder-chart.png', '_static/neowise-finder-chart.png')

If you choose a different rendering engine
instead of the venerable but rather ornery and complicated `matplotlib`_,
then of course the plotting calls you make
will be completely different.
But the basic data loading and filtering will be the same,
so hopefully the script will still help get you started
in targeting a more modern plotting library.

Plotting satellite altitude during re-entry
===========================================

Here is the decreasing altitude of a satellite as its orbit decayed
and it re-entered the atmosphere above the Pacific Ocean:

.. image:: _static/goce-reentry.png

The code to produce the diagram using `matplotlib`_,
including custom tick marks that are based on the date,
is:

.. testcode::

    from matplotlib import pyplot as plt
    from matplotlib.dates import HourLocator, DateFormatter

    from numpy import arange

    def label_dates_and_hours(axes):
        axes.xaxis.set_major_locator(HourLocator([0]))
        axes.xaxis.set_minor_locator(HourLocator([0, 12]))
        axes.xaxis.set_major_formatter(DateFormatter('\n%a %d'))
        axes.xaxis.set_minor_formatter(DateFormatter('%Hh'))

    from skyfield.api import load, EarthSatellite

    # Load the satellite's final TLE entry.

    sat = EarthSatellite(
        '1 34602U 09013A   13314.96046236  .14220718  20669-5  50412-4 0   930',
        '2 34602 096.5717 344.5256 0009826 296.2811 064.0942 16.58673376272979',
        'GOCE',
    )

    # Build the time range `t` over which to plot, plus other values.

    ts = load.timescale()
    t = ts.tt_jd(arange(sat.epoch.tt - 1.0, sat.epoch.tt + 3.0, 0.01))
    reentry = ts.utc(2013, 11, 11, 0, 16)
    earth_radius_km = 6371.

    # Start a new figure.

    fig, ax = plt.subplots()

    # Draw the blue curve.

    x = t.toordinal()
    y = sat.at(t).distance().km - earth_radius_km
    ax.plot(x, y)

    # Label the official moment of reentry.

    x = reentry.toordinal()
    y = sat.at(reentry).distance().km - earth_radius_km
    ax.plot(x, y, 'ro')
    ax.text(x, y + 10, 'Moment of re-entry')

    # Grid lines and labels.

    label_dates_and_hours(ax)
    ax.grid()
    ax.set(title='GOCE satellite altitude', ylabel='km above sea level')

    # Render the plot to a PNG file.

    fig.savefig('goce-reentry.png', bbox_inches='tight')

.. testcleanup::

    import os
    os.rename('goce-reentry.png', '_static/goce-reentry.png')

.. _matplotlib: https://matplotlib.org/
