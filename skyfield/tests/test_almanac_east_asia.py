import skyfield.almanac_east_asia
from skyfield import api, almanac

# Compare with Hong Kong Observatory:
# https://www.hko.gov.hk/tc/gts/astronomy/Solar_Term.htm access at 2019-12-14

def test_solar_terms():
    ts = api.load.timescale()
    e = api.load('de421.bsp')
    f = skyfield.almanac_east_asia.solar_terms(e)

    # https://en.wikipedia.org/wiki/Lichun

    t0 = ts.utc(2019, 2, 2)
    t1 = ts.utc(2019, 2, 5)
    t, y = almanac.find_discrete(t0, t1, f)
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2019-02-04 03:14']
    assert (y == (21)).all()

    # https://en.wikipedia.org/wiki/Lixia

    t0 = ts.utc(2019, 5, 4)
    t1 = ts.utc(2019, 5, 6)
    t, y = almanac.find_discrete(t0, t1, f)
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2019-05-05 19:03']
    assert (y == (3)).all()

    # https://en.wikipedia.org/wiki/Liqiu

    t0 = ts.utc(2019, 8, 6)
    t1 = ts.utc(2019, 8, 8)
    t, y = almanac.find_discrete(t0, t1, f)
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2019-08-07 19:13']
    assert (y == (9)).all()

    # https://en.wikipedia.org/wiki/Lidong

    t0 = ts.utc(2019, 11, 7)
    t1 = ts.utc(2019, 11, 9)
    t, y = almanac.find_discrete(t0, t1, f)
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2019-11-07 17:24']
    assert (y == (15)).all()
