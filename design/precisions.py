#!/usr/bin/env python

from skyfield.api import load
from math import tau

mas = tau / 360 / 3600 / 1e3
moon_distance_m = 384.4e6
venus_closest_m = 39.54e9
mars_closest_m = 62e9
mars_farthest_m = 401e9
jupiter_closest_m = 365e9
saturn_closest_m = 746e9

def main1():
    # See ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/de405.iom.pdf

    ts = load.timescale()
    de406 = load('de406.bsp')
    de422 = load('de422.bsp')
    de431 = load('de431t.bsp')

    t = ts.utc(-2990, range(1, 1 + 12 * (2990 * 2)))
    print(t[0].utc)
    print(t[-1].utc)

    #target = 'moon'
    target = 'jupiter barycenter'

    earth_moon_de406 = (de406[target] - de406['earth']).at(t)
    earth_moon_de422 = (de422[target] - de422['earth']).at(t)
    earth_moon_de431 = (de431[target] - de431['earth']).at(t)

    difference_406 = earth_moon_de431.separation_from(earth_moon_de406)
    difference_422 = earth_moon_de431.separation_from(earth_moon_de422)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(t.J, difference_422.arcseconds(), ',')
    ax.plot(t.J, difference_406.arcseconds(), ',')
    ax.grid()
    fig.savefig('tmp.png')

def main2():
    # "The error in the Earth and Mars orbits in DE 405 is now known to
    # be about 2 km" (Folkner, Williams, Boggs 2009)
    print('DE405, Mars:', to_mas(2e3 / mars_closest_m), 'mas')

    print('DE405, outer planets:', max(
        # Morrison and Evans 1998
        # "It is concluded that the overall accuracy of DE405 is better
        # than 50 mas for these planets in the period 1984–1997.
        50,
        # Stone 1998
        # "small offset, ∼0."09, however, in the right ascensions of Uranus
        90,
        # https://eclipse.gsfc.nasa.gov/SEhelp/limb.html
        # "the Jet Propusion Laboratory's solar system ephemeris (known
        # as DE405) is capable of Sun and Moon positions with errors of
        # less than 0.1 arc-seconds."
        100,
    ), 'mas')

    # Standish 1998:
    print('DE406 vs DE405:', max(
        # "no worse than 25 meters for all planets" + closest planet: venus
        to_mas(25.0 / venus_closest_m),
        # "no worse than 1 meter for the moon"
        to_mas(1.0 / moon_distance_m),
    ), 'mas')

    # Right half of residuals plot (Folkner, Williams, Boggs 2009)
    print('DE421, Moon:', to_mas(0.3 / moon_distance_m), 'mas')

    # Residuals plot (Folkner, Williams, Boggs 2009)
    print('DE421, Mercury:', to_mas(300 / venus_closest_m), 'mas')

    # "Venus orbit accuracy is now about 200 m" (Folkner, Williams, Boggs 2009)
    print('DE421, Venus:', to_mas(200 / venus_closest_m), 'mas')

    # "The Earth and Mars orbit accuracies are expected to be better
    # than 300 m through 2008" (Folkner, Williams, Boggs 2009)
    print('DE421, Mars:', to_mas(300 / mars_closest_m), 'mas')

    # "The orbits of Jupiter and Saturn are determined to accuracies of
    # tens of kilometers" plus residuals plot (Folkner, Williams, Boggs 2009)
    print('DE421, Jupiter:', to_mas(50e3 / jupiter_closest_m), 'mas')
    print('DE421, Saturn:', to_mas(50e3 / saturn_closest_m), 'mas')
    print('DE421, Saturn:', 500, 'mas')
    print('DE421, Uranus:', 1e3, 'mas')

    # (Folkner, Williams, Boggs, Park, and Kuchynka 2014)
    print('DE431 Moon:', max(
        # "present-day lunar orbit is known to submeter accuracy"
        to_mas(1.0 / moon_distance_m),
        0,
    ), 'mas')
    print()
    print('"Orbits of the inner planets are known to subkilometer accuracy"')
    print('DE431 0.0002 accuracy limit when Venus closest:',
          venus_closest_m * 0.2 * mas, 'm')
    print('DE431 0.0002 accuracy limit when Mars farthest:',
          mars_farthest_m * 0.2 * mas, 'm')

    print()
    print('DE431 outer planets:', max(
        #
        to_mas(1.0 / jupiter_closest_m),
        0,
    ), 'mas')

def to_mas(radians):
    return radians / tau * 360 * 3600 * 1000

"""

"Figure 1 shows Mars Odyssey range residuals relative to DE 418, which
was fit to range data through the end of 2006. DE 418 is seen to predict
the range to Mars 1 year into the future with an accuracy of about 15 m.
Similarly, DE 421 is expected to predict the Earth–Mars range to about
15 m through the end of 2008.  (The error in the plane-of-sky position
of Mars relative to Earth through the end of 2008 is about 300 m.)"

- https://ipnpr.jpl.nasa.gov/progress_report/42-178/178C.pdf

Principles, deduced from above and from elsewhere:

1. Error grows with time: they can only provide these bounds
   for a limited number of months in future.
   (When was paper written?  Dated August 15, 2009.)

2. Magnitude of residuals is NOT the same as the bound on future accuracy.

3. Error is not symmetric:
   we often know ranges with much greater accuracy
   than we know plane-of-sky positions,
   with the difference in this case running around a factor of 300m / 15m = 20.

"""

if __name__ == '__main__':
    main2()
