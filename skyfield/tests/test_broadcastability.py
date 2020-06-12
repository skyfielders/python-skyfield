import numpy as np
from skyfield.positionlib import ICRF

_deep = np.array([[[1]],[[0]],[[0]]])

def test_ecliptic_xyz_with_no_epoch():
    p = ICRF(_deep)
    x, y, z = p.ecliptic_xyz().au
    assert x.shape == y.shape == z.shape == (1, 1)
