import os
from skyfield.api import load, load_file

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

def test_multiple_non_overlapping_segments_per_target():
    ts = load.timescale()
    t = ts.utc(1969, [1, 4, 8, 12])
    eph = load_file(_data_path('./de441-1969.bsp'))
    pluto = eph['pluto barycenter']
    assert str(pluto.at(t).xyz) == (
        '[[-30.53621462 -30.51057842 -30.47219562 -30.42961491]\n'
        ' [ -0.31275432  -0.59574835  -0.97928763  -1.36269213]\n'
        ' [  9.10213562   9.00611101   8.87487322   8.74241232]] au'
    )

