from skyfield.api import Topos, load
from skyfield.constants import ERAD
from skyfield.framelib import itrs, true_equator_and_equinox_of_date
from skyfield.positionlib import Geocentric

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

def test_from_frame_method():
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)
    g1 = Geocentric([1,2,3], [4,5,6], t=t)
    r, v = g1.frame_xyz_and_velocity(itrs)  # which we trust: see the test above

    g2 = Geocentric.from_time_and_frame_vectors(t, itrs, r, v)
    assert max(abs(g2.position.au - [1,2,3])) < 2e-14
    assert max(abs(g2.velocity.au_per_d - [4,5,6])) < 3e-14

    # Make sure original vectors were not harmed (for example, by "+=").
    assert list(g1.position.au) == [1,2,3]
    assert list(g1.velocity.au_per_d) == [4,5,6]

def test_frame_without_spin():
    ts = load.timescale()
    t = ts.utc(2020, 11, 27, 15, 34)
    g = Geocentric([1,2,3], [4,5,6], t=t)

    # Simply test whether "None" spin raises an exception in either direction.
    f = true_equator_and_equinox_of_date
    r, v = g.frame_xyz_and_velocity(f)
    Geocentric.from_time_and_frame_vectors(t, f, r, v)
