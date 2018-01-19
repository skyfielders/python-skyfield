from skyfield.api import load
from skyfield.positionlib import Geocentric
from skyfield.toposlib import Topos

def ts():
    yield load.timescale()

def test_beneath(ts):
    t = ts.utc(2018, 1, 19, 14, 37, 55)
    def f(xyz): return str(Topos.beneath(Geocentric(xyz, None, t)))
    assert f([1, 0, 0]) == 'Topos 00deg 00\' 00.0" N 21deg 46\' 59.4" E'
    assert f([0, 0, 1]) == 'Topos 90deg 00\' 00.0" N 21deg 46\' 59.4" E'
    assert f([1, 0, 1]) == 'Topos 45deg 00\' 00.0" N 21deg 46\' 59.4" E'
