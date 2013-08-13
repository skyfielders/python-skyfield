"""Run a series of timing tests on core Skyfield atoms."""

from itertools import product
from numpy import array, einsum, mean, std, float64, zeros
from unittest import TestCase

from skyfield import (coordinates, earthlib, framelib, nutationlib,
                      planets, precessionlib, starlib, timescales)

from .constants import T0, DEG2RAD, AU_KM, TAU
from .timescales import julian_date, JulianDate
from timeit import default_timer

TA = julian_date(1969, 7, 20, 20. + 18. / 60.)
TB = julian_date(2012, 12, 21)

D0 = 63.8285
DA = 39.707
DB = 66.8779

import de421
emp = planets.Ephemeris(de421)
earth = emp.earth
star = starlib.Star(
    ra=1.59132070233, dec=8.5958876464,
    pm_ra=0.0, pm_dec=0.0,
    parallax=0.0, radial_velocity=0.0,
)


class BM(object):
    def __init__(self, times, bm_fn, t):
        self.times = times
        self.bm_fn = bm_fn
        self.t = t

    def __call__(self):
        self.bm_fn(self.times, self.t)


def run_benchmark(times, fn, *args, **kwargs):
    data = zeros(times)
    for i in xrange(times):
        start = default_timer()
        fn(*args, **kwargs)
        end = default_timer()

        data[i] = end - start

    avg, stdev = mean(data), std(data)
    suite_name = "{}.{}".format(fn.__module__, fn.__name__)
    print('{} : times = {}, avg = {}, std = {}'.format(
        suite_name, times, avg, stdev))


def bm_earth_rotation_angle(times, t):
    run_benchmark(times, earthlib.earth_rotation_angle, t)


def bm_star_observe_from(times, t):
    run_benchmark(times, star.observe_from, earth(t))


def bm_planet_observe_from(times, t):
    run_benchmark(times, emp.jupiter.observe_from, earth(t))


def bm_topo_planet_observe(times, t):
    ggr = coordinates.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.earth = emp.earth
    ggr.ephemeris = emp
    run_benchmark(times, ggr(t).observe, emp.jupiter)


def bm_earth_tilt(times, t):
    run_benchmark(times, nutationlib.earth_tilt, t)


def bm_equation_of_the_equinoxes(times, t):
    run_benchmark(times, nutationlib.equation_of_the_equinoxes_complimentary_terms, t)
    

def bm_coordinate_to_astrometric(times, t):
    coordinate = star.observe_from(earth(t))
    run_benchmark(times, coordinate.astrometric)


def bm_coordinate_to_apparent(times, t):
    coordinate = star.observe_from(earth(t))
    run_benchmark(times, coordinate.apparent)


def bm_coordinate_horizontal(times, t):
    ggr = coordinates.Topos('75 W', '45 N', 0.0,
                            temperature=10.0, pressure=1010.0)
    ggr.earth = emp.earth
    ggr.ephemeris = emp
    run_benchmark(times, ggr(t).observe(emp.jupiter).apparent().horizontal)


BENCHMARKS = (
    BM(times=100, bm_fn=bm_earth_rotation_angle, t=array([T0, TA, TB])),

    BM(times=100, bm_fn=bm_star_observe_from, t=JulianDate(tt=T0)),
    BM(times=100, bm_fn=bm_star_observe_from, t=JulianDate(tt=TA)),
    BM(times=100, bm_fn=bm_star_observe_from, t=JulianDate(tt=TB)),

    BM(times=100, bm_fn=bm_planet_observe_from, t=JulianDate(tt=T0)),
    BM(times=100, bm_fn=bm_planet_observe_from, t=JulianDate(tt=TA)),
    BM(times=100, bm_fn=bm_planet_observe_from, t=JulianDate(tt=TB)),

    BM(times=100, bm_fn=bm_topo_planet_observe, t=JulianDate(tt=T0)),
    BM(times=100, bm_fn=bm_topo_planet_observe, t=JulianDate(tt=TA)),
    BM(times=100, bm_fn=bm_topo_planet_observe, t=JulianDate(tt=TB)),

    BM(times=100, bm_fn=bm_earth_tilt, t=JulianDate(tt=T0)),
    BM(times=100, bm_fn=bm_earth_tilt, t=JulianDate(tt=TA)),
    BM(times=100, bm_fn=bm_earth_tilt, t=JulianDate(tt=TB)),

    BM(times=100, bm_fn=bm_equation_of_the_equinoxes, t=array([T0, TA, TB])),

    BM(times=100, bm_fn=bm_coordinate_to_astrometric, t=JulianDate(tt=T0)),
    BM(times=100, bm_fn=bm_coordinate_to_astrometric, t=JulianDate(tt=TA)),
    BM(times=100, bm_fn=bm_coordinate_to_astrometric, t=JulianDate(tt=TB)),

    BM(times=100, bm_fn=bm_coordinate_to_apparent, t=JulianDate(tt=T0)),
    BM(times=100, bm_fn=bm_coordinate_to_apparent, t=JulianDate(tt=TA)),
    BM(times=100, bm_fn=bm_coordinate_to_apparent, t=JulianDate(tt=TB)),

    BM(times=100, bm_fn=bm_coordinate_horizontal, t=JulianDate(tt=TB)),
)

if __name__ == "__main__":
    for bm in BENCHMARKS:
        bm()
