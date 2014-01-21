
==================
 Earth Satellites
==================

.. testcode::

    text = """
    ISS (ZARYA)             
    1 25544U 98067A   14020.93268519  .00009878  00000-0  18200-3 0  5082
    2 25544  51.6498 109.4756 0003572  55.9686 274.8005 15.49815350868473
    """

    from skyfield.api import earth

    bluffton = earth.topos('40.8939 N', '83.8917 W')
    #tup = (2014, 1, 21, 11, 18, 7)
    tup = (2014, 1, 21, 11, 18, 7)

    sat = earth.satellite(text)
    alt, az, distance = bluffton(utc=tup).observe(sat).altaz()

    print alt
    print az
    print distance.km

.. testoutput::

    13deg 31' 26.0"
    21deg 34' 12.0"
    1296.58733908


