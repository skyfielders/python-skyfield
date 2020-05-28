from numpy import array, isnan
from skyfield.geometry import intersect_line_and_sphere

def a(*args):
    return array(args)

def test_intersect_line_and_sphere():
    near, far = intersect_line_and_sphere(a(1000, 0, 0), a(2, 0, 0), 1)
    assert near == 1.0
    assert far == 3.0

    near, far = intersect_line_and_sphere(a(-1000, 0, 0), a(3, 0, 0), 1)
    assert near == -4.0
    assert far == -2.0

    near, far = intersect_line_and_sphere(a(0, 0.5, 0), a(0, 5, 0), 2)
    assert near == 3.0
    assert far == 7.0

    near, far = intersect_line_and_sphere(a(1000, 0, 0), a(0, 5, 0), 2)
    assert isnan(near)
    assert isnan(far)

    near, far = intersect_line_and_sphere(a(1000, 0, 0), a(2, 1, 0), 1)
    assert near == 2.0
    assert far == 2.0

    near, far = intersect_line_and_sphere(a(1000, 0, 0), a(2, 0.999, 0), 1)
    assert abs(near - 1.955) < 1e-3
    assert abs(far - 2.045) < 1e-3
