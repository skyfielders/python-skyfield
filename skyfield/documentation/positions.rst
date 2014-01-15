
===========================
 Positions and Coordinates
===========================

.. parsed-literal::

    Three positions

    obj(time)           →  Barycentric position (BCRS)
     └─ observe(obj2)   →  Astrometric position (ICRS)
         └─ apparent()  →  Apparent position (GCRS)

    Barycentric, Astrometric, or Apparent position
     │
     ├── `position <api.html#Position.position>`_.AU         →   x, y, z
     ├── `position <api.html#Position.position>`_.km         →   x, y, z
     ├── `position.to(unit) <api.html#Distance.to>`_   →   x, y, z
     │
     ├── `velocity <api.html#Position.velocity>`_.AU_per_d   →   xdot, ydot, zdot
     ├── `velocity <api.html#Position.Velocity>`_.km_per_s   →   xdot, ydot, zdot
     ├── `velocity.to(unit) <api.html#Distance.to>`_   →   xdot, ydot, zdot
     │
     ├── `radec() <api.html#Position.radec>`_             →   ra, dec, distance
     └── `radec(epoch=jd) <api.html#Position.radec>`_     →   ra, dec, distance

    Apparent position only
     │
     └── `altaz() <api.html#Position.altaz>`_             →   alt, az, distance

    Angle like ra, dec, alt, and az
     │
     ├── `radians() <api.html#Angle.radians>`_           →   6.266029488577352
     │
     ├── `hours() <api.html#Angle.hours>`_             →   23.934469599999996
     ├── `hms() <api.html#Angle.hms>`_               →   (1, 23, 56, 4, 0)
     ├── `hstr() <api.html#Angle.hstr>`_              →   '23h 56m 4.09s'
     ├── `hstr(places=4) <api.html#Angle.hstr>`_      →   '23h 56m 4.0906s'
     │
     ├── `degrees() <api.html#Angle.degrees>`_           →   359.017044
     ├── `dms() <api.html#Angle.dms>`_               →   (1, 359, 1, 1, 0)
     ├── `dstr() <api.html#Angle.dstr>`_              →   '359deg 1\' 1.4"'
     └── `dstr(places=3) <api.html#Angle.dstr>`_      →   '359deg 1\' 1.358"'

