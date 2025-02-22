"""Tests of how well we parse various file formats."""

import gzip
from skyfield import iokit
from skyfield.data import hipparcos, stellarium
from skyfield.functions import BytesIO
from skyfield.iokit import parse_tle

old_deltat_preds = b"""\
YEAR    TT-UT PREDICTION  UT1-UTC PREDICTION  ERROR

 2017.00      68.591             -0.408         0.000
 2017.25      68.72               0.469         0.00
 2017.50      68.81               0.376         0.01
 2017.75      68.86               0.322         0.01
 2018.00      68.99               0.192         0.02
 2018.25      69.14               0.041         0.02
 2018.50      69.3                              0.2
"""

def test_old_deltat_preds():
    lines = old_deltat_preds.splitlines()
    data = iokit.parse_deltat_preds(lines)

    assert data[0][0] == 2457754.5
    assert str(data[1][0]) == '68.591'

    assert data[0][-1] == 2458300.5
    assert str(data[1][-1]) == '69.3'

new_deltat_preds = b"""\
   MJD        YEAR    TT-UT Pred  UT1-UTC Pred  ERROR
   58484.000  2019.00   69.34      -0.152       0.117
   58575.000  2019.25   69.48      -0.295       0.162
   58666.000  2019.50   69.62      -0.440       0.215
   58758.000  2019.75   69.71      -0.527       0.273
   58849.000  2020.00   69.87                   0.335
"""

def test_new_deltat_preds():
    lines = new_deltat_preds.splitlines()
    data = iokit.parse_deltat_preds(lines)

    assert data[0][0] == 2458484.5
    assert str(data[1][0]) == '69.34'

    assert data[0][-1] == 2458849.5
    assert str(data[1][-1]) == '69.87'

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
    assert len(s) == 2
    assert s[0][0] == ['ISS (ZARYA)', 'ISS', 'ZARYA']
    assert s[0][1].name == 'ISS (ZARYA)'
    assert s[1][0] == ['FLOCK 2E-1']
    assert s[1][1].name == 'FLOCK 2E-1'

def test_spacetrack_two_line():
    f = BytesIO(sample_spacetrack_two_line_text)
    s = list(parse_tle(f))
    assert len(s) == 2
    assert s[0][0] == ()
    assert s[0][1].name is None
    assert s[1][0] == ()
    assert s[1][1].name is None

def test_spacetrack_three_line():
    f = BytesIO(sample_spacetrack_three_line_text)
    s = list(parse_tle(f))
    assert len(s) == 2
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
        df = hipparcos.load_dataframe(b)
    except ImportError:
        # raise SkipTest('pandas not available')
        # Assay doesn't understand skipping tests yet; just pass
        # for now if Pandas cannot be imported.
        return
    assert len(df) == 1
    row = df.iloc[0]
    assert abs(row.ra_degrees - 000.00091185) < 1e-30
    assert abs(row.dec_degrees - +01.08901332) < 1e-30

star_text = b"""\
# star names by constellation
# Andromeda (And)
   677|_("Alpheratz") 1,2,5,6,11,12
   677|_("Sirrah")
  5447|_("Mirach") 1,2,5,6,11,12,23
  9640|_("Almach") 1,2,5,6,11,12
  9640|_("Almaak")
"""

def test_stellarium_star_names():
    f = BytesIO(star_text)
    star_names = stellarium.parse_star_names(f)
    assert star_names[0].hip == 677
    assert star_names[0].name == 'Alpheratz'
    assert star_names[4].hip == 9640
    assert star_names[4].name == 'Almaak'
