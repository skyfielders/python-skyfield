from skyfield.api import Topos, load
from skyfield.constants import ERAD
from skyfield.framelib import itrs

def test_frame_rotation():
    # Does a frame's rotation and twist get applied in the right
    # directions?  Let's test whether the position and velocity of an
    # ITRS vector (ERAD,0,0) are restored to the proper orientation.
    top = Topos(latitude_degrees=0, longitude_degrees=0)
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)  # Arbitrary time; LST ~= 20.03.
    p = top.at(t)

    r = p.frame_xyz(itrs)
    assert max(abs(r.m - [ERAD, 0, 0])) < 4e-8 # meters

    r, v = p.frame_xyz_and_velocity(itrs)
    assert max(abs(r.m - [ERAD, 0, 0])) < 4e-8 # meters
    assert max(abs(v.km_per_s)) < 3e-15 # km/s
