from numpy import add, append, array, diff, flatnonzero, reshape, sign

def _remove_adjacent_duplicates(indices):
    mask = diff(indices) != 0
    mask = append(mask, [True])
    return indices[mask]

def _choose_brackets(y):
    dsd = diff(sign(diff(y)))
    indices = flatnonzero(dsd < 0)
    left = reshape(add.outer(indices, [0, 1]), -1)
    left = _remove_adjacent_duplicates(left)
    right = left + 1
    return left, right

def test_brackets_of_simple_peak():
    y = array((10, 11, 12, 11, 10))
    assert [list(v) for v in _choose_brackets(y)] == [[1, 2], [2, 3]]

def test_brackets_of_small_plateau():
    y = array((10, 11, 12, 12, 11, 10))
    assert [list(v) for v in _choose_brackets(y)] == [[1, 2, 3], [2, 3, 4]]

def test_brackets_of_wide_plateau():
    y = array((10, 11, 12, 12, 12, 12, 12, 11, 10))
    assert [list(v) for v in _choose_brackets(y)] == [[1, 2, 5, 6],
                                                      [2, 3, 6, 7]]
