"""Determine whether arrays work as well as individual inputs."""

import sys
from numpy import array
from .. import starlib
from ..constants import T0, B1950
from ..planets import earth, mars
from ..positionlib import Topos
from ..timescales import JulianDate, julian_date

if sys.version_info < (3,):
    from itertools import izip
else:
    izip = zip

dates = array([
    julian_date(1969, 7, 20, 20. + 18. / 60.),
    T0,
    julian_date(2012, 12, 21),
    julian_date(2027, 8, 2, 10. + 7. / 60. + 50. / 3600.),
    ])

deltas = array([39.707, 63.8285, 66.8779, 72.])

def observe_planet_from_geocenter(ut1, delta_t):
    jd = JulianDate(ut1=ut1, delta_t=delta_t)

    yield jd.ut1
    yield jd.tt
    yield jd.tdb

    observer = earth(jd)

    yield observer.position
    yield observer.velocity
    yield observer.jd.ut1
    yield observer.jd.tt
    yield observer.jd.tdb

    astrometric = observer.observe(mars)

    yield astrometric.position
    yield astrometric.velocity

    ra, dec, distance = astrometric.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance

    ra, dec, distance = astrometric.radec(epoch=B1950)

    yield ra.hours()
    #yield dec.degrees()
    #yield distance

    apparent = astrometric.apparent()

    yield apparent.position
    #yield apparent.velocity  # = None?

    ra, dec, distance = apparent.radec()

    yield ra.hours()
    #yield dec.degrees()
    #yield distance

    ra, dec, distance = apparent.radec(epoch=B1950)

    yield ra.hours()
    #yield dec.degrees()
    #yield distance

def observe_planet_from_topos(ut1, delta_t):
    jd = JulianDate(ut1=ut1, delta_t=delta_t)

    yield jd.ut1
    yield jd.tt
    yield jd.tdb

    topos = Topos('71.1375 W', '42.6583 N', 0.0)
    topos.ephemeris = earth.ephemeris
    observer = topos(jd)

    yield observer.position
    yield observer.velocity
    yield observer.jd.ut1
    yield observer.jd.tt
    yield observer.jd.tdb

    astrometric = observer.observe(mars)

    yield astrometric.position
    yield astrometric.velocity

    ra, dec, distance = astrometric.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance

    ra, dec, distance = astrometric.radec(epoch=B1950)

    yield ra.hours()
    #yield dec.degrees()
    #yield distance

    apparent = astrometric.apparent()

    yield apparent.position
    #yield apparent.velocity  # = None?

    ra, dec, distance = apparent.radec()

    yield ra.hours()
    #yield dec.degrees()
    #yield distance

    ra, dec, distance = apparent.radec(epoch=B1950)

    yield ra.hours()
    #yield dec.degrees()
    #yield distance

def compute_stellar_position(ut1, delta_t):
    star = starlib.Star(ra=1.59132070233, dec=8.5958876464)
    observer = earth(ut1=ut1, delta_t=delta_t)
    astrometric = observer.observe(star)

    yield astrometric.position
    yield astrometric.velocity

    ra, dec, distance = astrometric.radec()

    yield ra.hours()
    yield dec.degrees()
    yield distance

def pytest_generate_tests(metafunc):
    if 'vector_vs_scalar' in metafunc.fixturenames:
        metafunc.parametrize(
            'vector_vs_scalar',
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
    assert vector is not None, (
        '{}:\n  vector is None'.format(location))
    assert (vector.T[i] == scalar).all(), (
        '{}:\n  {}[{}] != {}\n  difference: {}'.format(
            location, vector.T, i, scalar, vector.T[i] - scalar))
