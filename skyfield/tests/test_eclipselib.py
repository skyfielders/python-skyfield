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
