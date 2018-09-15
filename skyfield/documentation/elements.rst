
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
:class:`~skyfield.elementslib.OsculatingElements`.  For example, here is
how to find the osculating elements of the moon orbiting earth:

 .. testcode::

    from skyfield.api import load
    from skyfield.elementslib import osculating_elements_of

    ts = load.timescale()
    t = ts.utc(2018, 4, 22, range(0, 25))

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

Attributes of OsculatingElements objects
========================================

Here is a list of the attributes of the Elements object and their types:

.. parsed-literal::

    Element describing shape of the orbit:
     └── ``eccentricity``                → numpy.ndarray

    Element describing the tilt of the orbital plane:
     └── ``inclination``                 → Angle object

    Element describing the direction in which the orbital plane is tilted:
     └── ``longitude_of_ascending node`` → Angle object

    Elements describing direction of periapsis:
     ├── ``argument_of_periapsis``       → Angle object
     ├── ``longitude_of_periapsis``      → Angle object
     └── ``periapsis_time``              → Time object

    Elements describing the size of the orbit:
     ├── ``apoapsis_distance``           → Distance object
     ├── ``mean_motion_per_day``         → Angle object
     ├── ``periapsis_distance``          → Distance object
     ├── ``period_in_days``              → numpy.ndarray
     ├── ``semi_latus_rectum``           → Distance object
     ├── ``semi_major_axis``             → Distance object
     └── ``semi_minor_axis``             → Distance object

    Elements describing the secondary's position in the orbit:
     ├── ``argument_of_latitude``        → Angle object
     ├── ``eccentric_anomaly``           → Angle object
     ├── ``mean_anomaly``                → Angle object
     ├── ``mean_longitude``              → Angle object
     ├── ``true_anomaly``                → Angle object
     ├── ``true_longitude``              → Angle object
     └── (the secondary's position can be implicit in periapsis_time 
              because at periapsis all anomalies are 0)

    Other attributes:
     └── ``time``                        → Time object

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
