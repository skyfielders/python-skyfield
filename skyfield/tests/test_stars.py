import numpy as np
from assay import assert_raises
from skyfield import api
from skyfield.data.hipparcos import load_dataframe
from skyfield.errors import EphemerisRangeError

def eph():
    eph = api.load('de421.bsp')
    try:
        yield eph
    finally:
        eph.close()

def test_dataframe():
    with api.load.open('hip_main.dat.gz') as f:
        df = load_dataframe(f)
    star = api.Star.from_dataframe(df)
    assert repr(star) == 'Star(ra shape=9933, dec shape=9933, ra_mas_per_year shape=9933, dec_mas_per_year shape=9933, parallax_mas shape=9933, epoch shape=9933)'

def test_apparent_nan_star(eph):
    ts = api.load.timescale()
    earth = eph['earth']
    surface = earth + api.wgs84.latlon(35.0, -111.0)
    for name in ('ra_hours', 'dec_degrees', 'ra_mas_per_year',
                 'dec_mas_per_year', 'parallax_mas', 'radial_km_per_s'):
        kwargs = dict(ra_hours=1.0, dec_degrees=2.0)
        kwargs[name] = np.nan
        star = api.Star(**kwargs)
        for t in (ts.utc(2026, 5, 19), ts.utc(2026, 5, [19, 20])):
            for observer in (earth, surface):
                with np.errstate(all='raise'):
                    apparent = observer.at(t).observe(star).apparent()
                    ra, dec, distance = apparent.radec()
                assert apparent.xyz.au.shape == (3,) + t.shape
                assert np.isnan(apparent.xyz.au).all()
                assert np.isnan(ra.hours).all()
                assert np.isnan(dec.degrees).all()
                assert np.isnan(distance.au).all()

def test_apparent_mixed_nan_stars(eph):
    ts = api.load.timescale()
    t = ts.utc(2026, 5, 19)
    earth = eph['earth']
    surface = earth + api.wgs84.latlon(35.0, -111.0)
    ra = np.array([1.0, np.nan, 2.0, 3.0])
    dec = np.array([4.0, 5.0, 6.0, np.nan])
    good = np.array([True, False, True, False])
    stars = api.Star(ra_hours=ra, dec_degrees=dec)
    valid_stars = api.Star(ra_hours=ra[good], dec_degrees=dec[good])

    for observer in (earth, surface):
        observer = observer.at(t)
        astrometric = observer.observe(stars)
        original = astrometric.xyz.au.copy()
        with np.errstate(all='raise'):
            apparent = astrometric.apparent()
            expected = observer.observe(valid_stars).apparent()
        assert apparent.xyz.au.shape == (3, 4)
        assert np.isnan(apparent.xyz.au[:, ~good]).all()
        np.testing.assert_array_equal(apparent.xyz.au[:, good], expected.xyz.au)
        np.testing.assert_array_equal(astrometric.xyz.au, original)
        if observer._observer_gcrs_au is not None:
            with np.errstate(all='raise'):
                alt, az, distance = apparent.altaz()
            assert np.isnan(alt.degrees[~good]).all()
            assert np.isnan(az.degrees[~good]).all()
            assert np.isnan(distance.au[~good]).all()
            assert np.isfinite(alt.degrees[good]).all()
            assert np.isfinite(az.degrees[good]).all()

def test_apparent_all_nan_and_empty_star_arrays(eph):
    ts = api.load.timescale()
    observer = eph['earth'].at(ts.utc(2026, 5, 19))
    for n in (0, 1, 3):
        values = np.full(n, np.nan)
        stars = api.Star(ra_hours=values, dec_degrees=values)
        with np.errstate(all='raise'):
            apparent = observer.observe(stars).apparent()
        assert apparent.xyz.au.shape == (3, n)
        assert np.isnan(apparent.xyz.au).all()

def test_apparent_dataframe_with_nan_stars(eph):
    with api.load.open('hip_main.dat.gz') as f:
        df = load_dataframe(f)
    ts = api.load.timescale()
    observer = eph['earth'].at(ts.utc(2026, 5, 19))
    good = df['ra_degrees'].notnull() & df['dec_degrees'].notnull()
    assert (~good).sum() == 5
    with np.errstate(all='raise'):
        apparent = observer.observe(api.Star.from_dataframe(df)).apparent()
        expected = observer.observe(api.Star.from_dataframe(df[good])).apparent()
    assert apparent.xyz.au.shape == (3, len(df))
    assert np.isnan(apparent.xyz.au[:, ~good]).all()
    np.testing.assert_array_equal(apparent.xyz.au[:, good], expected.xyz.au)

def test_apparent_still_raises_for_out_of_range_deflector_time(eph):
    ts = api.load.timescale()
    # The observer is in range, but light passing the Sun before the
    # kernel's first instant still needs to raise a real range error.
    t = ts.tdb_jd(eph.segments[0].spk_segment.start_jd)
    observer = eph['earth'].at(t)
    sun = (eph['sun'] - eph['earth']).at(t)
    ra, dec, _ = sun.radec()
    star = api.Star(ra=ra, dec=dec)
    with assert_raises(EphemerisRangeError):
        observer.observe(star).apparent()
