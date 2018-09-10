from skyfield import api, almanac

# http://aa.usno.navy.mil/cgi-bin/aa_moonill2.pl?form=1&year=2018&task=00&tz=-05
def test_fraction_illuminated():
    ts = api.load.timescale()
    t0 = ts.utc(2018, 9, range(9, 19), 5)
    e = api.load('de421.bsp')
    p = almanac.fraction_illuminated(e, 'moon', t0).round(2)
    assert (p == (0, 0, 0.03, 0.08, 0.15, 0.24, 0.33, 0.43, 0.52, 0.62)).all()
