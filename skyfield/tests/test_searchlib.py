from numpy import array, int8
from skyfield.api import load
from skyfield.searchlib import (
    _choose_brackets, _fix_numpy_deprecation, _identify_maxima, find_discrete,
)

def test_brackets_of_simple_peak():
    y = array((10, 11, 12, 11, 10))
    left, right = _choose_brackets(y)
    assert list(left) == [1, 2]
    assert list(right) == [2, 3]

def test_brackets_of_small_plateau():
    y = array((10, 11, 12, 12, 11, 10))
    left, right = _choose_brackets(y)
    assert list(left) == [1, 2, 3]
    assert list(right) == [2, 3, 4]

def test_brackets_of_wide_plateau():
    y = array((10, 11, 12, 12, 12, 12, 12, 11, 10))
    left, right = _choose_brackets(y)
    assert list(left) == [1, 2, 5, 6]
    assert list(right) == [2, 3, 6, 7]

def test_simple_maxima():
    x = array((2451545.0, 2451546.0, 2451547.0))
    y = array((11, 12, 11))
    x, y = _identify_maxima(x, y)
    assert list(x) == [2451546.0]
    assert list(y) == [12]

def test_maxima_of_small_plateau():
    x = array((2451545.0, 2451546.0, 2451547.0, 2451548.0))
    y = array((11, 12, 12, 11))
    x, y = _identify_maxima(x, y)
    assert list(x) == [2451546.5]
    assert list(y) == [12]

def test_both_kinds_of_maxima_at_once():
    # We put what could be maxima at each end, to make sure they are not
    # detected as when we can't confirm them.
    y = array((12, 11, 12, 12, 11, 13, 11, 12, 12))
    x = 2451545.0 + array(range(len(y)))
    x, y = _identify_maxima(x, y)
    assert list(x) == [2451547.5, 2451550.0]
    assert list(y) == [12, 13]

def test_boolean_values_become_int8():
    y = _fix_numpy_deprecation(array((True, False, True)))
    assert y.dtype == int8
    assert list(y) == [1, 0, 1]

def test_find_discrete_with_boolean_function():
    # NumPy 2.5 warns if code sets an array's `dtype` (issue #1153).
    ts = load.timescale()
    def f(t):
        return t.tt % 1.0 > 0.5
    f.step_days = 0.1
    t, y = find_discrete(ts.tt_jd(2460000.0), ts.tt_jd(2460002.0), f)
    assert y.dtype == int8
    assert list(y) == [1, 0, 1, 0]
