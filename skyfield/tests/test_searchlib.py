from numpy import (add, append, argsort, array, concatenate, diff,
                   flatnonzero, reshape, sign)

def _remove_adjacent_duplicates(indices):
    mask = diff(indices) != 0
    mask = append(mask, [True])
    return indices[mask]

def _choose_brackets(y):
    """Return the indices between which we should search for maxima of `y`."""
    dsd = diff(sign(diff(y)))
    indices = flatnonzero(dsd < 0)
    left = reshape(add.outer(indices, [0, 1]), -1)
    left = _remove_adjacent_duplicates(left)
    right = left + 1
    return left, right

def _identify_maxima(x, y):
    dsd = diff(sign(diff(y)))

    # Choose every point that is higher than the two adjacent points.
    indices = flatnonzero(dsd == -2) + 1
    peak_x = x.take(indices)
    peak_y = y.take(indices)

    # Also choose the midpoint between the edges of a plateau, if both
    # edges are in view.  First we eliminate runs of zeroes, then look
    # for adjacent -1 values, then map those back to the main array.
    indices = flatnonzero(dsd)
    dsd2 = dsd.take(indices)
    minus_ones = dsd2 == -1
    plateau_indices = flatnonzero(minus_ones[:-1] & minus_ones[1:])
    plateau_left_indices = indices.take(plateau_indices)
    plateau_right_indices = indices.take(plateau_indices + 1) + 2
    plateau_x = x.take(plateau_left_indices) + x.take(plateau_right_indices)
    plateau_x /= 2.0
    plateau_y = y.take(plateau_left_indices + 1)

    x = concatenate((peak_x, plateau_x))
    y = concatenate((peak_y, plateau_y))
    indices = argsort(x)
    return x[indices], y[indices]

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
