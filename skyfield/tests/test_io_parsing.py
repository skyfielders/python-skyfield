"""Tests of how well we parse various file formats."""

from skyfield.functions import BytesIO
from skyfield.iokit import parse_celestrak_tle

sample_celestrak_text = b"""\
ISS (ZARYA)             \n\
1 25544U 98067A   18135.61844383  .00002728  00000-0  48567-4 0  9998
2 25544  51.6402 181.0633 0004018  88.8954  22.2246 15.54059185113452
FLOCK 2E-1              \n\
1 41483U 98067JD  18135.38689952  .00096183  14684-4  28212-3 0  9990
2 41483  51.6270 103.3896 0004826  61.7810 298.3684 15.92672255114129
"""

def test_celestrak():
    f = BytesIO(sample_celestrak_text)
    d = dict(parse_celestrak_tle(f))
    assert len(d) == 4
    assert d['ISS'] == d['ISS (ZARYA)'] == d['ZARYA']
    assert d['FLOCK 2E-1']
