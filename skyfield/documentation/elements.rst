
=============================
 Osculating Orbital Elements
=============================

.. currentmodule:: skyfield.api

Skyfield is able to calculate osculating orbital elements for any
Geometric, Barycentric, or Geocentric position with respect to either
the ecliptic or equatorial plane. The data produced by Skyfield matches
the data produced by JPL’s HORIZONS system.

Generating Elements
===================

Call :func:`~skyfield.elementslib.osculating_elements_of()` to generate
an :class:`~skyfield.elementslib.OsculatingElements` object.  For
example, here is how to find the osculating elements of the moon
orbiting earth:

.. testcode::

    from skyfield.api import load
    from skyfield.elementslib import osculating_elements_of

    ts = load.timescale()
    t = ts.utc(2018, 4, 22)

    planets = load('de421.bsp')
    earth = planets['earth']
    moon = planets['moon']

    position = (moon - earth).at(t)
    elements = osculating_elements_of(position)

The elements are then attributes of the Elements object:

.. testcode::

    i = elements.inclination.degrees
    e = elements.eccentricity
    a = elements.semi_major_axis.km

    print('Inclination: {0:.2f} degrees'.format(i))
    print('Eccentricity: {0:.5f}'.format(e))
    print('Semimajor axis: {0:.0f} km'.format(a))

.. testoutput::

    Inclination: 20.46 degrees
    Eccentricity: 0.03104
    Semimajor axis: 380577 km

Note that one element,
the periapsis time,
is not necessarily a unique value:
if an orbit is periodic,
then the body will reach periapsis repeatedly
at a series of dates which are separated by ``period_in_days``.

.. testcode::

    print('Periapsis:', elements.periapsis_time.utc_strftime())
    print('Period: {0:.2f} days'.format(elements.period_in_days))

.. testoutput::

    Periapsis: 2018-04-20 16:09:42 UTC
    Period: 26.88 days

You can add or subtract the period
in order to produce a series of equally valid periapsis dates
for that set of orbital elements.

.. testcode::

    next = elements.periapsis_time + elements.period_in_days
    print('Next periapsis:', next.utc_strftime())

.. testoutput::

    Next periapsis: 2018-05-17 13:14:56 UTC

Attributes of OsculatingElements objects
========================================

Here is a list of the attributes of the Elements object and their types:

.. parsed-literal::

    **OsculatingElements object**
     │   **Element describing the shape of the orbit:**
     ├── eccentricity                → `numpy.ndarray <https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html>`_
     │
     │   **Element describing the tilt of the orbital plane:**
     ├── inclination                 → `Angle object <api-units.html>`_
     │
     │   **Element describing the direction in which the orbital plane is tilted:**
     ├── longitude_of_ascending_node → `Angle object <api-units.html>`_
     │
     │   **Elements describing direction of periapsis:**
     ├── argument_of_periapsis       → `Angle object <api-units.html>`_
     ├── longitude_of_periapsis      → `Angle object <api-units.html>`_
     ├── periapsis_time              → `Time object <api.html#time-objects>`_
     │
     │   **Elements describing the size of the orbit:**
     ├── apoapsis_distance           → `Distance object <api-units.html>`_
     ├── mean_motion_per_day         → `Angle object <api-units.html>`_
     ├── periapsis_distance          → `Distance object <api-units.html>`_
     ├── period_in_days              → `numpy.ndarray <https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html>`_
     ├── semi_latus_rectum           → `Distance object <api-units.html>`_
     ├── semi_major_axis             → `Distance object <api-units.html>`_
     ├── semi_minor_axis             → `Distance object <api-units.html>`_
     │
     │   **Elements describing the secondary's position in the orbit:**
     ├── argument_of_latitude        → `Angle object <api-units.html>`_
     ├── eccentric_anomaly           → `Angle object <api-units.html>`_
     ├── mean_anomaly                → `Angle object <api-units.html>`_
     ├── mean_longitude              → `Angle object <api-units.html>`_
     ├── true_anomaly                → `Angle object <api-units.html>`_
     ├── true_longitude              → `Angle object <api-units.html>`_
     ├── (the secondary's position can be implicit in periapsis_time 
     │        because at periapsis all anomalies are 0)
     │
     │   **Other attributes:**
     └── time                        → `Time object <api.html#time-objects>`_

To fully define an object's location and orbit, one element is required from each of the above categories.

Reference Planes
================

By default the ``elements()`` method produces elements using the *xy*
plane of the ICRF as the reference plane.  This is equivalent to the
J2000.0 equatorial plane within the tolerance of J2000.0.  If you
instead want elements using the J2000.0 ecliptic as the reference plane,
pass it as the second argument:

 .. testcode::

    from skyfield.data.spice import inertial_frames
    ecliptic = inertial_frames['ECLIPJ2000']

    t = ts.utc(2018, 4, 22, range(0,25))
    position = (moon - earth).at(t)
    elements = osculating_elements_of(position, ecliptic)
