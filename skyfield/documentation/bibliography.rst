
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

    http://aa.usno.navy.mil/software/novas/novas_f/NOVAS_F3.1_Guide.pdf

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

*Delta T: Past, Present and Future* at the United States Naval Observatory

    http://asa.usno.navy.mil/SecK/DeltaT.html

    A few well-presented plots of ∆T and its behavior over the past few
    centuries.  I tended to have this page always open in another tab
    while reading about ∆T on other sites, so that I could correlate the
    descriptions of each text against the visuals here.

*IERS Rapid Service/Prediction Center*
links at the United States Naval Observatory

    http://maia.usno.navy.mil/

    The many links on this page helped me sort out the various raw
    sources of Earth orientation data that are available, and decide on
    the ones that needed to be built into Skyfield.  I also made many
    visits to the IERS web site itself:

    http://www.iers.org/IERS/EN/DataProducts/EarthOrientationData/eop.html

    but found the USNO site to always be the better starting point.

*The ∆T pages at the NASA Eclipse Web Site*

    | http://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html

    While most of the pages about ∆T at the Eclipse Web Site are a basic
    introduction to the concept, the “Polynomial Expressions” page
    provides something new: the polynomials that Espenak and Meeus used
    as their fit to ∆T when building their *Five Millennium Canon of
    Solar Eclipses: −1999 to +3000.*

    When the day comes that I implement eclipse logic in Skyfield, I
    will be interested in comparing my results against the `Espenak and
    Meeus paper itself`_!

.. _NOVAS library: http://aa.usno.navy.mil/software/novas/novas_py/novaspy_intro.php
.. _Espenak and Meeus paper itself: http://eclipse.gsfc.nasa.gov/5MCSE/5MCSE-Text11.pdf
.. _jplephem: https://pypi.python.org/pypi/jplephem
