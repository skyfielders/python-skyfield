"""Run a series of timing tests on core Skyfield atoms."""

from itertools import product
from numpy import array, einsum, mean, std, float64, zeros
from unittest import TestCase

from skyfield import (coordinates, earthlib, framelib, nutationlib,
                      planets, precessionlib, starlib, timescales)

from .constants import T0, DEG2RAD, AU_KM, TAU
from .timescales import julian_date
from timeit import default_timer


def benchmark(times, fn, *args, **kwargs):
    data = zeros(times)
    for i in xrange(times):
        start = default_timer()
        fn(*args, **kwargs)
        end = default_timer()

        data[i] = end - start

    avg, stdev = mean(data), std(data)
    print('{} : times = {}, avg = {}, std = {}'.format(
        fn.__module__, times, avg, stdev))


def test_earth_rotation_angle():
    TA = julian_date(1969, 7, 20, 20. + 18. / 60.)
    TB = julian_date(2012, 12, 21)
    t = array([T0, TA, TB])

    benchmark(1000, earthlib.earth_rotation_angle, t)

if __name__ == "__main__":
    test_earth_rotation_angle()
