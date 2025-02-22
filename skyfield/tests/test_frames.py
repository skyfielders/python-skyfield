from numpy import cos
from skyfield import framelib
from skyfield.api import SSB, Star, Topos, load, wgs84
from skyfield.constants import AU_M, ERAD
from skyfield.positionlib import Geocentric

def test_radec_and_altaz_angles_and_rates():
    # HORIZONS test data in Skyfield repository: authorities/radec-altaz-rates
    ts = load.timescale()
    t = ts.utc(2021, 2, 3)
    top = wgs84.latlon(35.1844866, 248.347300, elevation_m=2106.9128)
    planets = load('de421.bsp')
    a = (planets['earth'] + top).at(t).observe(planets['mars']).apparent()

    # First, verify RA and declination.

    frame = framelib.true_equator_and_equinox_of_date
    dec, ra, distance, dec_rate, ra_rate, range_rate = (
        a.frame_latlon_and_rates(frame))

    arcseconds = 3600.0
    assert abs((ra.degrees - 40.75836) * arcseconds) < 0.04
    assert abs((dec.degrees - 17.16791) * arcseconds) < 0.005
    assert abs(distance.m - 1.21164331503552 * AU_M) < 120.0

    # Verify RA and declination rates of change.

    assert round(dec_rate.arcseconds.per_hour, 5) == 25.61352
    assert round(ra_rate.arcseconds.per_hour * cos(dec.radians),
                 4) == round(75.15571, 4)  # TODO: get last digit to agree?
    assert abs(range_rate.km_per_s - 16.7926932) < 2e-5

    # Verify altitude and azimuth.

    frame = top
    alt, az, distance, alt_rate, az_rate, range_rate = (
        a.frame_latlon_and_rates(frame))

    assert round(alt.degrees, 4) == 65.2758
    assert round(az.degrees, 4) == 131.8839
    assert abs(distance.m - 1.21164331503552 * AU_M) < 120.0

    # Verify altitude and azimuth rates of change.

    assert abs(range_rate.km_per_s - 16.7926932) < 2e-5
    assert round(alt_rate.arcseconds.per_minute, 2) == 548.66
    assert round(az_rate.arcseconds.per_minute * cos(alt.radians), 2) == 663.55

def test_frame_round_trip():
    # Does a frame's rotation and twist get applied in the right
    # directions?  Let's test whether the position and velocity of an
    # ITRS vector (ERAD,0,0) are restored to the proper orientation.
    top = Topos(latitude_degrees=0, longitude_degrees=0)
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)  # Arbitrary time; LST ~= 20.03.
    p = top.at(t)

    r = p.frame_xyz(framelib.itrs)
    assert max(abs(r.m - [ERAD, 0, 0])) < 4e-8 # meters

    r, v = p.frame_xyz_and_velocity(framelib.itrs)
    assert max(abs(r.m - [ERAD, 0, 0])) < 4e-8 # meters
    assert max(abs(v.km_per_s)) < 3e-15 # km/s

def test_from_frame_method():
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)
    g1 = Geocentric([1,2,3], [4,5,6], t=t)
    r, v = g1.frame_xyz_and_velocity(framelib.itrs) # which we trust: see above

    g2 = Geocentric.from_time_and_frame_vectors(t, framelib.itrs, r, v)
    assert max(abs(g2.xyz.au - [1,2,3])) < 2e-14
    assert max(abs(g2.velocity.au_per_d - [4,5,6])) < 3e-14

    # Make sure original vectors were not harmed (for example, by "+=").
    assert list(g1.xyz.au) == [1,2,3]
    assert list(g1.velocity.au_per_d) == [4,5,6]

def test_frame_without_spin():
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)
    g = Geocentric([1,2,3], [4,5,6], t=t)

    # Simply test whether "None" spin raises an exception in either direction.
    f = framelib.true_equator_and_equinox_of_date
    r, v = g.frame_xyz_and_velocity(f)
    Geocentric.from_time_and_frame_vectors(t, f, r, v)

def test_true_equator_and_equinox_of_date():
    ts = load.timescale()
    t = ts.utc(2025, 7, 2.375)
    alpha_andromeda = Star(
        ra_hours=0.13976888866666667, dec_degrees=29.09082805,
        ra_mas_per_year=135.68, dec_mas_per_year=-162.95, parallax_mas=33.6,
        epoch=ts.J(1991.25),
    )
    dec, ra, _ = SSB.at(t).observe(alpha_andromeda).frame_latlon(
        framelib.mean_equator_and_equinox_of_date
    )
    ra.preference = 'hours'

    # Position from page H2 of the 2025 Astronomical Almanac:
    assert ra.hstr(places=1) == '00h 09m 42.7s'
    assert dec.dstr() == '+29deg 13\' 52.0"'

def test_tirs_at_least_runs():
    # TODO: find an external source for a TIRS vector to test against.
    # For now, just make sure it doesn't raise an exception.
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)
    g = Geocentric([1,2,3], [4,5,6], t=t)
    g.frame_xyz_and_velocity(framelib.tirs)
