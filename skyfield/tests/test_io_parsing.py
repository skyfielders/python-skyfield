"""Tests of how well we parse various file formats."""

from skyfield.functions import BytesIO
from skyfield.iokit import parse_tle

sample_celestrak_text = b"""\
ISS (ZARYA)             \n\
1 25544U 98067A   18135.61844383  .00002728  00000-0  48567-4 0  9998
2 25544  51.6402 181.0633 0004018  88.8954  22.2246 15.54059185113452
FLOCK 2E-1              \n\
1 41483U 98067JD  18135.38689952  .00096183  14684-4  28212-3 0  9990
2 41483  51.6270 103.3896 0004826  61.7810 298.3684 15.92672255114129
"""

sample_spacetrack_text = b"""\
1 29273U 06033B   18081.29838594 -.00000056 +00000-0 +00000-0 0  9993
2 29273 000.0189 154.5198 0004980 202.4902 284.9321 01.00271755042548
1 29274U 06033C   18081.39999693 +.00002637 +00000-0 +10299-2 0  9992
2 29274 005.9144 244.7152 6177908 248.3941 037.5897 03.74556424124616
"""

def test_celestrak():
    f = BytesIO(sample_celestrak_text)
    s = list(parse_tle(f))
    print(s)
    assert len(s) == 2
    assert s[0][0] == ['ISS (ZARYA)', 'ISS', 'ZARYA']
    assert s[0][1].name == 'ISS (ZARYA)'
    assert s[1][0] == ['FLOCK 2E-1']
    assert s[1][1].name == 'FLOCK 2E-1'

def test_spacetrack():
    f = BytesIO(sample_spacetrack_text)
    s = list(parse_tle(f))
    assert len(s) == 2
    print(s)
    assert s[0][0] == ()
    assert s[0][1].name is None
    assert s[1][0] == ()
    assert s[1][1].name is None

def test_extra_lines_in_tle_are_ignored():
    f = BytesIO(b'Sample line\n' + sample_celestrak_text + b'Another line\n')
    s = list(parse_tle(f))
    assert len(s) == 2
    assert s[0][0] == ['ISS (ZARYA)', 'ISS', 'ZARYA']
    assert s[0][1].name == 'ISS (ZARYA)'
    assert s[1][0] == ['FLOCK 2E-1']
    assert s[1][1].name == 'FLOCK 2E-1'
