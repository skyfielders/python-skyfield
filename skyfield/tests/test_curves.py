from skyfield.curvelib import Splines, build_spline_given_ends

def test_spline_and_derivative():
    parameters = 10, 12, 2.0, 3.0, 5.0, 7.0
    x = 8.0, 9.0, 10.0, 11.0, 12.0
    expected_values = 3.0, 5.0, 7.0, 10.5, 17.0
    expected_slope = 2.5, 1.75, 2.5, 4.75, 8.5

    # Does the Spline class work?

    curve = Splines(parameters)
    assert tuple(curve(x)) == expected_values
    assert tuple(curve.derivative(x)) == expected_slope

    # Can we rebuild the parameters from its endpoints and slopes?

    parameters2 = build_spline_given_ends(
        10, curve(10), curve.derivative(10),
        12, curve(12), curve.derivative(12))

    assert parameters == parameters2

    # If we rebuild the spline around different bounds, does it still
    # produce the same curve?

    for lower, upper in (8,10), (8,9), (21,23):
        parameters3 = build_spline_given_ends(
            lower, curve(lower), curve.derivative(lower),
            upper, curve(upper), curve.derivative(upper))
        curve3 = Splines(parameters3)
        assert tuple(curve3(x)) == expected_values
        assert tuple(curve3.derivative(x)) == expected_slope
