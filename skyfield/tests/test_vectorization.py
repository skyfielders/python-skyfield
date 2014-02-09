"""Determine whether arrays work as well as individual inputs."""

import sys
from numpy import array, rollaxis
from .. import starlib
from ..constants import T0, B1950
from ..api import earth, mars
from ..positionlib import Topos
from ..timelib import JulianDate, julian_date

if sys.version_info < (3,):
    from itertools import izip
else:
    izip = zip

dates = array([
    julian_date(1969, 7, 20) + (20.0 + 18.0 / 60.0) / 24.0,
    T0,
    julian_date(2012, 12, 21),
    julian_date(2027, 8, 2) + (10.0 + (7.0 + 50.0 / 60.0) / 60.0) / 24.0,
    ])

deltas = array([39.707, 63.8285, 66.8779, 72.])

def compute_times_and_equinox_matrices(tt, delta_t):
    jd = JulianDate(tt=tt, delta_t=delta_t)

    yield jd.ut1
    yield jd.tt
    yield jd.tdb

    yield jd.P
    yield jd.N
    yield jd.M

def observe_planet_from_geocenter(tt, delta_t):
    jd = JulianDate(tt=tt, delta_t=delta_t)
    observer = earth(jd)

    yield observer.position.AU
    yield observer.position.km
    yield observer.velocity.AU_per_d
    yield observer.velocity.km_per_s
    yield observer.jd.ut1
    yield observer.jd.tt
    yield observer.jd.tdb

    astrometric = observer.observe(mars)

    yield astrometric.position.AU
    yield astrometric.velocity.AU_per_d

    ra, dec, distance = astrometric.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

    ra, dec, distance = astrometric.radec(epoch=B1950)

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

    apparent = astrometric.apparent()

    yield apparent.position.AU
    #yield apparent.velocity  # = None?

    ra, dec, distance = apparent.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

    ra, dec, distance = apparent.radec(epoch=B1950)

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

def observe_planet_from_topos(tt, delta_t):
    jd = JulianDate(tt=tt, delta_t=delta_t)

    yield jd.ut1
    yield jd.tt
    yield jd.tdb

    topos = Topos('71.1375 W', '42.6583 N', 0.0)
    topos.ephemeris = earth.ephemeris
    observer = topos(jd)

    yield observer.position.AU
    yield observer.velocity.AU_per_d
    yield observer.jd.ut1
    yield observer.jd.tt
    yield observer.jd.tdb

    astrometric = observer.observe(mars)

    yield astrometric.position.AU
    yield astrometric.velocity.AU_per_d

    ra, dec, distance = astrometric.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

    ra, dec, distance = astrometric.radec(epoch=B1950)

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

    apparent = astrometric.apparent()

    yield apparent.position.AU
    #yield apparent.velocity  # = None?

    ra, dec, distance = apparent.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

    ra, dec, distance = apparent.radec(epoch=B1950)

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

def compute_stellar_position(tt, delta_t):
    star = starlib.Star(ra=1.59132070233, dec=8.5958876464)
    observer = earth(tt=tt, delta_t=delta_t)
    astrometric = observer.observe(star)

    yield astrometric.position.AU
    yield astrometric.velocity.AU_per_d

    ra, dec, distance = astrometric.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance.AU

def pytest_generate_tests(metafunc):
    if 'vector_vs_scalar' in metafunc.fixturenames:
        metafunc.parametrize(
            'vector_vs_scalar',
            list(generate_comparisons(compute_times_and_equinox_matrices)) +
            list(generate_comparisons(observe_planet_from_geocenter)) +
            list(generate_comparisons(observe_planet_from_topos)) +
            list(generate_comparisons(compute_stellar_position)))

def generate_comparisons(computation):
    """Set up comparisons between vector and scalar outputs of `computation`.

    The `computation` should be a generator that accepts both vector and
    scalar input, and that yields a series of values whose shape
    corresponds to its input's shape.

    """
    vector_results = list(computation(dates, deltas))
    for i, (date, delta_t) in enumerate(zip(dates, deltas)):
        g = computation(date, delta_t)
        for vector, scalar in izip(vector_results, g):
            f = g.gi_frame
            location = '{}:{}'.format(f.f_code.co_filename, f.f_lineno)
            yield location, vector, i, scalar

def test_vector_vs_scalar(vector_vs_scalar):
    location, vector, i, scalar = vector_vs_scalar
    vectorT = rollaxis(vector, -1)
    assert vector is not None, (
        '{}:\n  vector is None'.format(location))
    assert vectorT[i].shape == scalar.shape, (
        '{}:\n  {}[{}].shape != {}.shape\n  shapes: {} {}'.format(
            location, vector.T, i, scalar, vector.T[i].shape, scalar.shape))

    vectorTi = vectorT[i]

    # Yes, an auto-generated epsilon with no physical significance!
    # Why?  Because we are comparing the rounding differences in two
    # (hopefully!) identical floating-point computations, not thinking
    # of the results as two physical calculations.

    epsilon = abs(1e-15 * max(vectorTi.max(), scalar.max()))
    difference = abs(vectorTi - scalar)

    assert (difference <= epsilon).all(), (
        '{}:\n vector[{}] = {}\n'
        ' scalar    = {}\n'
        ' difference= {}\n'
        ' epsilon   = {}'
        .format(location, i, vectorTi, scalar, difference, epsilon))
