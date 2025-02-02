# Some tests of the behavior of the almanac `_intersection()` function.

def _q(a, b, c, sign):
    # when a x^2 + b x + c = 0
    from numpy import sqrt
    # print('doing quadratic with:', a, b, c)
    print('---')
    print('root (+):', (-b + sqrt(b*b - 4*a*c)) / 2*a)
    print('root (-):', (-b - sqrt(b*b - 4*a*c)) / 2*a)
    print('root alt (+):', - 2*c / (b + sqrt(b*b - 4*a*c)))
    print('root alt (-):', - 2*c / (b - sqrt(b*b - 4*a*c)))
    print('sign:', sign)
    return - 2*c / (b + sign * sqrt(b*b - 4*a*c))

def _intersection(a0, a1, v0, v1):
    # Return the time at which a curve reaches a=0, given its position
    # and velocity a0, v0 at time 0.0 and a1, v1 at time 1.0.
    #
    # (overdetermined, so, ignores v1)
    # print('intersection with:', a0, a1, v0, v1)
    # print('k would be:', 2 * (a1 - a0 - v0))
    sign = 1 - 2 * (a0 > a1)
    tx = _q(a1 - a0 - v0, v0, a0, sign)
    # print('tx =', tx)
    return tx

def test_intersection():
    for v in 0.9, 1, 1.1:
        print(v)
        assert _intersection(0, 4, v, v) == 0.0
        assert _intersection(-5, 0, v, v) == 1.0

    assert _intersection(-1, 3, 4, 4) == 0.25
    assert _intersection(-3, 1, 4, 4) == 0.75

    # Increasing, with other intersection far to the right.

    t = _intersection(-5, 5, 10.1, 9.9)
    assert t > 0.49
    assert t < 0.50

    # Increasing, with other intersection far to the left.

    t = _intersection(-5, 5, 9.9, 10.1)
    assert t > 0.50
    assert t < 0.51

    # Decreasing, with other intersection far to the right.

    t = _intersection(5, -5, -10.1, -9.9)
    assert t > 0.49
    assert t < 0.50

    # Decreasing, with other intersection far to the left.

    t = _intersection(5, -5, -9.9, -10.1)
    assert t > 0.50
    assert t < 0.51
    asdf
