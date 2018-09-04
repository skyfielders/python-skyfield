"""Tests of how well we parse various file formats."""

import gzip
from skyfield.data.hipparcos import load_dataframe
from skyfield.functions import BytesIO
from skyfield.iokit import parse_tle

sample_celestrak_text = b"""\
ISS (ZARYA)             \n\
1 25544U 98067A   18135.61844383  .00002728  00000-0  48567-4 0  9998
2 25544  51.6402 181.0633 0004018  88.8954  22.2246 15.54059185113452
FLOCK 2E-1              \r\n\
1 41483U 98067JD  18135.38689952  .00096183  14684-4  28212-3 0  9990
2 41483  51.6270 103.3896 0004826  61.7810 298.3684 15.92672255114129
"""

sample_spacetrack_two_line_text = b"""\
1 29273U 06033B   18081.29838594 -.00000056 +00000-0 +00000-0 0  9993
2 29273 000.0189 154.5198 0004980 202.4902 284.9321 01.00271755042548
1 29274U 06033C   18081.39999693 +.00002637 +00000-0 +10299-2 0  9992
2 29274 005.9144 244.7152 6177908 248.3941 037.5897 03.74556424124616
"""

sample_spacetrack_three_line_text = b"""\
0 First
1 29273U 06033B   18081.29838594 -.00000056 +00000-0 +00000-0 0  9993
2 29273 000.0189 154.5198 0004980 202.4902 284.9321 01.00271755042548
0 Second
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

def test_spacetrack_two_line():
    f = BytesIO(sample_spacetrack_two_line_text)
    s = list(parse_tle(f))
    assert len(s) == 2
    print(s)
    assert s[0][0] == ()
    assert s[0][1].name is None
    assert s[1][0] == ()
    assert s[1][1].name is None

def test_spacetrack_three_line():
    f = BytesIO(sample_spacetrack_three_line_text)
    s = list(parse_tle(f))
    assert len(s) == 2
    print(s)
    assert s[0][0] == ['First']
    assert s[0][1].name == 'First'
    assert s[1][0] == ['Second']
    assert s[1][1].name == 'Second'

def test_extra_lines_in_tle_are_ignored():
    f = BytesIO(b'Sample line\n' + sample_celestrak_text + b'Another line\n')
    s = list(parse_tle(f))
    assert len(s) == 2
    assert s[0][0] == ['ISS (ZARYA)', 'ISS', 'ZARYA']
    assert s[0][1].name == 'ISS (ZARYA)'
    assert s[1][0] == ['FLOCK 2E-1']
    assert s[1][1].name == 'FLOCK 2E-1'

sample_hipparcos_line = b"""\
H|           1| |00 00 00.22|+01 05 20.4| 9.10| |H|000.00091185|+01.08901332| |   3.54|   -5.20|   -1.88|  1.32|  0.74|  1.39|  1.36|  0.81| 0.32|-0.07|-0.11|-0.24| 0.09|-0.01| 0.10|-0.01| 0.01| 0.34|  0| 0.74|     1| 9.643|0.020| 9.130|0.019| | 0.482|0.025|T|0.55|0.03|L| | 9.2043|0.0020|0.017| 87| | 9.17| 9.24|       | | | |          | |  | 1| | | |  |   |       |     |     |    |S| | |224700|B+00 5077 |          |          |0.66|F5          |S \n\
"""

def test_hipparcos():
    b = BytesIO()
    g = gzip.GzipFile(mode='wb', fileobj=b)
    g.write(sample_hipparcos_line)
    g.close()
    b.seek(0)
    try:
        df = load_dataframe(b)
    except ImportError:
        # raise SkipTest('pandas not available')
        # Assay doesn't understand skipping tests yet; just pass
        # for now if Pandas cannot be imported.
        return
    assert len(df) == 1
    row = df.iloc[0]
    assert abs(row.ra_degrees - 000.00091185) < 1e-30
    assert abs(row.dec_degrees - +01.08901332) < 1e-30
