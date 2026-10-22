import numpy as np
from skyfield import api


def ts():
    yield api.load.timescale()


def test_indexed_observation_can_be_made_apparent(ts):
    eph = api.load('de421.bsp')
    observer = eph['earth'] + api.wgs84.latlon(37.7749, -122.4194)
    times = ts.utc(2025, 1, [1, 2, 3])
    position = observer.at(times).observe(eph['mars'])
    expected = position.apparent()
    for index in [1, slice(1, None), [2, 0, 2], np.array([True, False, True])]:
        actual = position[index].apparent()
        np.testing.assert_array_equal(actual.xyz.au, expected.xyz.au[:, index])
        np.testing.assert_array_equal(position[index].light_time,
                                      position.light_time[index])
        np.testing.assert_array_equal(position[index].phase_angle(eph['sun']).radians,
                                      position.phase_angle(eph['sun']).radians[index])
        np.testing.assert_array_equal(actual.altaz()[0].radians,
                                      expected.altaz()[0].radians[index])


def test_indexed_observer_can_observe(ts):
    eph = api.load('de421.bsp')
    observer = eph['earth'] + api.wgs84.latlon(37.7749, -122.4194)
    times = ts.utc(2025, 1, [1, 2, 3])
    position = observer.at(times)
    expected = position.observe(eph['mars']).apparent()
    for index in [1, slice(1, None), [2, 0, 2]]:
        selected = position[index]
        assert selected.center_barycentric is selected
        actual = selected.observe(eph['mars']).apparent()
        np.testing.assert_array_equal(actual.xyz.au, expected.xyz.au[:, index])


def test_indexed_star_array_keeps_scalar_observer(ts):
    eph = api.load('de421.bsp')
    observer = (eph['earth'] + api.wgs84.latlon(37.7749, -122.4194)).at(ts.utc(2025))
    stars = api.Star(ra_hours=[1., 2., 3.], dec_degrees=[10., 20., 30.])
    position = observer.observe(stars)
    expected = position.apparent()
    for index in [1, slice(1, None), [2, 0, 2]]:
        selected = position[index]
        assert selected.center_barycentric is observer
        np.testing.assert_array_equal(selected.apparent().xyz.au,
                                      expected.xyz.au[:, index])


def test_indexed_satellite_position_keeps_earth_occultation(ts):
    eph = api.load('de421.bsp')
    satellite = api.EarthSatellite(
        '1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993',
        '2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106')
    times = ts.utc(2018, 7, 4, 0, range(10))
    observer = eph['earth'] + satellite
    position = observer.at(times).observe(eph['moon']).apparent()
    expected = position.is_behind_earth()
    for index in [1, slice(1, None), [2, 0, 2]]:
        np.testing.assert_array_equal(position[index].is_behind_earth(), expected[index])
