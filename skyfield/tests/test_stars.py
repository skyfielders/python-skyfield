from skyfield import api
from skyfield.data.hipparcos import load_dataframe
from skyfield.framelib import ICRS

def test_astrometric_velocity():
    ts = api.load.timescale()
    for zero, jd in ((0.0, api.T0),
                     (0.0, [api.T0, api.T0 + 1.0]),
                     ([0.0, 0.0], api.T0)):
        star = api.Star(ra_hours=zero, dec_degrees=zero,
                        ra_mas_per_year=1000.0, parallax_mas=1000.0)
        a = api.SSB.at(ts.tdb_jd(jd)).observe(star)
        ra_rate = a.frame_latlon_and_rates(ICRS)[4]
        assert (abs(ra_rate.mas.per_day - 1000.0 / 365.25) < 1e-9).all()

def test_dataframe():
    with api.load.open('hip_main.dat.gz') as f:
        df = load_dataframe(f)
    star = api.Star.from_dataframe(df)
    assert repr(star) == 'Star(ra shape=9933, dec shape=9933, ra_mas_per_year shape=9933, dec_mas_per_year shape=9933, parallax_mas shape=9933, epoch shape=9933)'
