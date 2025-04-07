import os
from skyfield.api import load, load_file

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
