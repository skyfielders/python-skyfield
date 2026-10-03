from skyfield.api import load
from skyfield import eclipselib

def test_lunar_eclipses():
    # The documentation test already confirms the dates of these two
    # eclipses; here, we confirm that the data structures all match.
    ts = load.timescale()
    eph = load('de421.bsp')

    t0 = ts.utc(2019, 1, 1)
    t1 = ts.utc(2020, 1, 1)
    t, y, details = eclipselib.lunar_eclipses(t0, t1, eph)

    assert len(t) == len(y) == 2
    for name, item in details.items():
        assert len(item) == len(t)

def test_solar_eclipses():
    ts = load.timescale()
    eph = load('de421.bsp')

    t0 = ts.utc(2019, 1, 1)
    t1 = ts.utc(2020, 1, 1)
    times, codes, mag, obs = eclipselib.solar_eclipses(t0, t1, eph)

    result = []
    for t, c, m, o in zip(times, codes, mag, obs):
        result += [f"{t.tt_strftime()} {eclipselib.SOLAR_ECLIPSES[c]} {m:.4} {o:.2}"]

    assert ','.join(result) == ("2019-01-06 01:42:40 TT Partial 0.7192 0.63,"
                                "2019-07-02 19:24:04 TT Total/Hybrid 1.046 1.0,"
                                "2019-12-26 05:18:52 TT Annular 0.9699 0.94")
