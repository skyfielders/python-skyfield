import os
from skyfield.api import load, load_file

from assay import assert_raises

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

# Test file generated with:
# python -m jplephem excerpt 1969/07/29 1969/07/30 de441.bsp de441-1969.bsp

def test_multiple_non_overlapping_segments_per_target():
    ts = load.timescale()
    t = ts.utc(1969, 7, [28, 29, 30, 31])
    eph = load_file(_data_path('./de441-1969.bsp'))

    pluto = eph['pluto barycenter']
    assert type(pluto).__name__ == 'Stack'
    assert len(pluto.segments) == 2

    assert str(pluto.at(t).xyz) == (
        '[[-30.47352052 -30.47318972 -30.47285863 -30.47252727]\n'
        ' [ -0.96671442  -0.96985773  -0.97300104  -0.97614434]\n'
        ' [  8.87919555   8.87811509   8.87703455   8.87595392]] au'
    )

    # Does the Stack correctly offer its .ephemeris to apparent()?

    pluto.at(t).observe(pluto).apparent()

    # TODO: SSB.at(t).observe() fails the above test.

# Verify that ephemeris objects let their segments be edited.

def test_removing_segments_from_jpl_ephemeris():
    eph = load('de421.bsp')
    eph.segments = [s for s in eph.segments if s.target in (3, 4)]

    assert len(eph.segments) == 2
    assert 2 not in eph
    assert 3 in eph
    assert eph.codes == {0, 3, 4}
    assert eph.names() == {
        0: ['SOLAR_SYSTEM_BARYCENTER', 'SSB', 'SOLAR SYSTEM BARYCENTER'],
        3: ['EARTH_BARYCENTER', 'EMB', 'EARTH MOON BARYCENTER',
            'EARTH-MOON BARYCENTER', 'EARTH BARYCENTER'],
        4: ['MARS_BARYCENTER', 'MARS BARYCENTER'],
    }

    assert repr(eph) == "<SpiceKernel 'de421.bsp'>"
    assert str(eph) == """\
Segments from kernel file 'de421.bsp':
  JD 2414864.50 - JD 2471184.50  (1899-07-28 through 2053-10-08)
      0 -> 3    SOLAR SYSTEM BARYCENTER -> EARTH BARYCENTER
      0 -> 4    SOLAR SYSTEM BARYCENTER -> MARS BARYCENTER"""

    assert type(eph[4]).__name__ == 'ChebyshevPosition'
    with assert_raises(KeyError):
        eph[5]
