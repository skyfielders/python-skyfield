
===============
 Example Plots
===============

This section of the documentation
will gradually accumulate example scripts
for producing images from Skyfield computations.
For the moment there’s only example so far,
for plotting the elevation of a satellite over time:

Finder chart for comet NEOWISE
==============================

.. testsetup::

    import matplotlib
    matplotlib.use('Agg')  # to avoid “no display name” error on Travis CI
    del matplotlib

.. testcode::

    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt

    from skyfield import api
    from skyfield.api import load
    ts = load.timescale(builtin=True)
    eph = load('de421.bsp')
    sun = eph['sun']
    earth = eph['earth']


    # In[3]:


    t_comet = ts.utc(2020, 7, range(17, 27))
    t = t_comet[len(t_comet) // 2]  # middle date


    # In[4]:


    from skyfield.data import mpc

    with load.open(mpc.COMET_URL) as f:
        comets = mpc.load_comets_dataframe(f)

    comets = comets.set_index('designation', drop=False)
    row = comets.loc['C/2020 F3 (NEOWISE)']


    # In[5]:


    from skyfield.constants import GM_SUN_Pitjeva_2005_km3_s2 as GM_SUN

    comet = sun + mpc.comet_orbit(row, ts, GM_SUN)
    center = earth.at(t).observe(comet)


    # In[6]:


    from skyfield.api import Star
    from skyfield.data import hipparcos


    # In[53]:


    from skyfield.projections import build_stereographic_projection
    proj = build_stereographic_projection(center)


    # In[54]:


    with load.open(hipparcos.URL) as f:
        stars = hipparcos.load_dataframe(f)


    # In[55]:


    star_positions = earth.at(t).observe(Star.from_dataframe(stars))
    stars['x'], stars['y'] = proj(star_positions)


    # In[56]:


    limiting_magnitude = 6.5
    bright_stars = (stars.magnitude <= limiting_magnitude)


    # In[57]:


    comet_x, comet_y = proj(earth.at(t_comet).observe(comet))


    # In[98]:


    # = 'https://raw.githubusercontent.com/Stellarium/stellarium/master/skycultures/western_SnT/constellationship.fab'
    from skyfield.data.stellarium import parse_constellations

    url = 'https://raw.githubusercontent.com/Stellarium/stellarium/master/skycultures/western_SnT/constellationship.fab'
    with load.open(url) as f:
        constellations = parse_constellations(f)

    edges_star1 = [star1 for name, edges in constellations for star1, star2 in edges]
    edges_star2 = [star2 for name, edges in constellations for star1, star2 in edges]

    np.array([stars['x'].loc[edges_star1], stars['y'].loc[edges_star1]])
    xy1 = stars[['x', 'y']].loc[edges_star1].values
    xy2 = stars[['x', 'y']].loc[edges_star2].values
    lines_xy = np.rollaxis(np.array([xy1, xy2]), 1)


    # In[94]:


    fig, ax = plt.subplots(figsize=[9, 9])

    from matplotlib.collections import LineCollection
    from skyfield.api import tau

    field_of_view_degrees = 45.0

    line_collection = LineCollection(lines_xy) #, color=['#0002'] * len(lines_xy))
    ax.add_collection(line_collection)

    marker_size = (0.5 + limiting_magnitude - stars.magnitude[bright_stars]) ** 2.0
    ax.scatter(stars['x'][bright_stars], stars['y'][bright_stars], s=marker_size, color='k')
    #ax.scatter(comet_x, comet_y, s=100, color='b')
    #comet_color = '#f00'
    #ax.plot(comet_x, comet_y, '+', c=comet_color, zorder=3)
    ax.plot(comet_x, comet_y, '+', zorder=3)
    offset = 0.002
    for xi, yi, tstr in zip(comet_x, comet_y, t_comet.utc_strftime('%-m/%d')):
        #text = ax.text(xi + offset, yi - offset, tstr, color=comet_color, ha='left', va='top', fontsize=12,
        text = ax.text(xi + offset, yi - offset, tstr, ha='left', va='top', fontsize=12,
               weight='bold')
        text.set_alpha(0.3)

    ax.set_title('Comet NEOWISE {} through {}'.format(
        t_comet[0].utc_strftime('%Y %B %d'),
        t_comet[-1].utc_strftime('%Y %B %d'),
    ))
    ax.set_aspect(1.0)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    angle = (180.0 - field_of_view_degrees) / 720.0 * tau
    limit = np.sin(angle) / (1.0 - np.cos(angle))
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    fig.savefig('neowise-finder-chart.png')



.. image:: _static/neowise-finder-chart.png

.. testcleanup::

    import os
    os.rename('neowise-finder-chart.png', '_static/neowise-finder-chart.png')

Satellite altitude during re-entry
==================================

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
