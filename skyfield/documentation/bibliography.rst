
==============
 Bibliography
==============

Further reading on the concepts and technology behind Skyfield.

*USNO Circular 179*
by George H. Kaplan, October 2005

    https://www.usno.navy.mil/USNO/astronomical-applications/publications/Circular_179.pdf/at_download/file

    My favorite introduction to how the International Astronomical Union
    responded in 1997 and 2000 to the problems with the old J2000 system
    by adopting the ICRS and a new Earth rotation model.  Thanks to
    Kaplan’s explanations, I finally understood the maneuvers of the
    USNO’s freely available `NOVAS library`_ (which comes with a Python
    interface!)  when dealing with dates and coordinate systems, and I
    was able to perform the same calculations correctly in Skyfield.  It
    also provided a useful introduction to the JPL ephemerides that
    NOVAS uses.

*User’s Guide to NOVAS Version F3.1,* March 2011

    https://aa.usno.navy.mil/downloads/novas/NOVAS_C3.1_Guide.pdf

    There were also several useful tidbits of information in the United
    States Naval Observatory’s guide to using NOVAS.  It tends towards
    more practical and user oriented advice than does the Circular, and
    of course was also my guide to the library’s concepts when I used it
    as my model for Skyfield’s internal astrometry routines.

*SPICE Required Reading documents*

    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/kernel.html
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
    https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html

    When implementing Skyfield’s sister library `jplephem`_, these
    documents were invaluable for their detailed description of the .bsp
    binary file format as well as the description of how the SPICE
    system combines multiple ephemeris segments to learn the
    displacement between two Solar System objects.

*JPL Solar System Dynamics directory of ephemeris PDFs*

    ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/

    I cannot always find the official documentation describing a
    particular JPL ephemeris, but this directory is always a good
    starting point, as it collects several of the PDFs together in one
    place.

.. _Morrison, Stephenson, et al:

*Earth Rotation — the Change in the Length of Day
and ΔT Plot showing lod from -2000 to +2500*

    http://astro.ukho.gov.uk/nao/lvm/

    This page offers the most up-to-date ΔT tables
    from Morrison, Stephenson, Hohenkerk, and Zawilski,
    which Skyfield uses to predict the Earth’s orientation
    for the years preceding and following the more detailed
    numbers published for 1973 through the current day by the IERS (see below).
    It also links to their academic papers,
    which provide a wealth of information
    about the historical events — primarily eclipses —
    that provide us with evidence
    for which direction the Earth was pointing in past centuries.

*Historical values of the Earth's clock error ∆T and the calculation of
eclipses* by Morrison & Stephenson, 2004

    | http://adsabs.harvard.edu/full/2004JHA....35..327M
    | `Full PDF <http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?2004JHA....35..327M&data_type=PDF_HIGH&type=PRINTER&filetype=.pdf>`_

    Detailed discussion of the fact that ∆T has the shape of a parabola
    over history because, omitting short-term variations, the shortening
    length of Earth’s day makes the difference between modern clock time
    and actual sunrise and sunset skew ever longer over the centuries as
    the tiny error introduced each year gradually accumulates.  Ancient
    eclipses for which records survive are our one point-source of data
    about how far the error had accumulated each century.

*IERS “Earth orientation data” page*

    There used to be a United States Naval Observatory web page
    dedicated to Earth orientation,
    but these days I look for data files here:

    http://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html

.. _NOVAS library: https://aa.usno.navy.mil/software/novaspy_intro
.. _jplephem: https://pypi.python.org/pypi/jplephem
