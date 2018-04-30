
=============================
 Osculating Orbital Elements
=============================

.. currentmodule:: skyfield.api

Skyfield is able to calculate osculating orbital elements for any 
Geometric, Barycentric, or Geocentric position with respect to either 
the ecliptic or equatorial plane. The data produced by skyfield matches 
the data produced by JPL's Horizons system.

Generating Elements
===================

An ``OsculatingElements`` object is returned by the ``elements()`` method of Geometric, 
Barycentric, and Geocentric positions. For example, here is how to find the 
osculating elements of the moon orbiting earth:

  .. testcode::
  
    t = ts.utc(2018, 4, 22, arange(0,25))
    elements = (moon - earth).at(t).elements()
    
The elements are then attributes of the Elements object:

  .. testcode::
  
    i = elem.inclination.degrees
    e = elem.eccentricity
    a = elem.semi_major_axis.km

Attributes of OsculatingElements objects
========================================

Here is a list of the attributes of the Elements object and their types:

.. parsed-literal::

    Elements describing orientation of the orbit plane:
     ├── ``inclination``                 → Angle object
     └── ``longitude_of_periapsis``      → Angle object
     
    Elements describing shape and direction the orbit within its plane: 
     ├── ``argument_of_periapsis``       → Angle object
     ├── ``eccentricity``                → ndarray
     ├── ``longitude_of_periapsis``      → Angle object
     └── ``periapsis_time``              → Time object
     
    Elements describing the size of the orbit:
     ├── ``apoapsis_distance``           → Distance object
     ├── ``mean_motion_per_day``         → ndarray
     ├── ``periapsis_distance``          → Distance object
     ├── ``period_in_days``              → ndarray
     ├── ``semi_latus_rectum``           → Distance object
     ├── ``semi_major_axis``             → Distance object
     └── ``semi_minor_axis``             → Distance object
     
    Elements describing the secondary's position in the orbit:
     ├── ``argument_of_latitude``        → Angle object
     ├── ``eccentric_anomaly``           → Angle object
     ├── ``mean_anomaly``                → Angle object
     ├── ``mean_longitude``              → Angle object
     ├── ``true_anomaly``                → Angle object
     └── ``true_longitude``              → Angle object
     
    Other attributes:
     └── ``time``                        → Time object
     
Reference Planes
================

By default the ``elements()`` method produces elements using 
the XY plane of the ICRF as the reference plane. This is 
equivalent to the J2000.0 equatorial plane within the tolerance 
of J2000.0. If you instead want elements using the J2000.0 ecliptic 
as the reference plane set the keyword argument ``ref_plane`` 
to ``'ecliptic'`` like this:

  .. testcode::
  
    t = ts.utc(2018, 4, 22, arange(0,25))
    elements = (moon - earth).at(t).elements(ref_plane='ecliptic')
