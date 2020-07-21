
===============
 Example Plots
===============

This section of the documentation
will gradually accumulate example scripts
for producing images from Skyfield computations.
For the moment there’s only example so far,
for plotting the elevation of a satellite over time:

.. testsetup::

    import matplotlib
    matplotlib.use('Agg')  # to avoid “no display name” error on Travis CI
    del matplotlib

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

    ts = load.timescale(builtin=True)
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

    fig.savefig('goce-reentry.png')

.. image:: _static/goce-reentry.png

.. testcleanup::

    import os
    os.rename('goce-reentry.png', '_static/goce-reentry.png')
