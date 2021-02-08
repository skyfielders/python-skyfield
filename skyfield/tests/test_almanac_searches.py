"""Low-level tests of the almanac search routines."""

import numpy as np

from skyfield.api import load
from skyfield.searchlib import find_discrete, find_maxima

bump = 1e-5
epsilon = 1e-10

def make_t():
    ts = load.timescale()
    t0 = ts.tt_jd(0)
    t1 = ts.tt_jd(1)
    return t0, t1

def make_stairstep_f(steps):
    """Return a function that increases by one at each of several `steps`."""
    def f(t):
        # For each time, sum how many of the values in `steps` it surpasses.
        return np.greater_equal.outer(t.tt, steps).sum(axis=1)
    f.step_days = 0.3
    return f

def is_close(value, expected):
    return (abs(value - expected) < epsilon).all()

def test_find_discrete_that_finds_nothing():
    t0, t1 = make_t()
    f = make_stairstep_f([-0.1, +1.1])
    t, y = find_discrete(t0, t1, f, epsilon)
    assert not len(t.tt)
    assert not len(y)

def test_find_discrete_near_left_edge():
    t0, t1 = make_t()
    f = make_stairstep_f([bump, 0.5])
    t, y = find_discrete(t0, t1, f, epsilon)
    assert is_close(t.tt, (bump, 0.5))
    assert list(y) == [1, 2]

def test_find_discrete_near_right_edge():
    t0, t1 = make_t()
    f = make_stairstep_f([0.5, 1.0 - bump])
    t, y = find_discrete(t0, t1, f, epsilon)
    assert is_close(t.tt, (0.5, 1.0 - bump))
    assert list(y) == [1, 2]

def test_find_discrete_with_a_barely_detectable_jag_right_at_zero():
    t0, t1 = make_t()
    f = make_stairstep_f([0.5, 0.5 + 3.1 * epsilon])
    t, y = find_discrete(t0, t1, f, epsilon)
    assert is_close(t.tt, (0.5, 0.5 + 3.1 * epsilon))
    assert list(y) == [1, 2]

def test_find_discrete_with_a_sub_epsilon_jag_right_at_zero():
    t0, t1 = make_t()
    f = make_stairstep_f([0.5, 0.5 + 0.99 * epsilon])

    # We hard-code num=12, just in case the default ever changes to
    # another value that might not trigger the symptom.
    t, y = find_discrete(t0, t1, f, epsilon, 12)

    # Note that we always return the last of several close solutions, so
    # that `y` correctly reflects the new state that persists after the
    # flurry of changes is complete.
    assert is_close(t.tt, (0.5 + 0.99 * epsilon,))
    assert list(y) == [2]

def make_mountain_range_f(peaks):
    """Return a function with local maxima at each of a series of `peaks`."""
    def f(t):
        # For each time, sum how many of the values in `steps` it surpasses.
        return -abs(np.subtract.outer(t.tt, peaks)).min(axis=1)
    f.step_days = 0.3
    return f

def test_finding_enough_maxima():
    # If the step size is small enough, no maxima should be skipped.
    t0, t1 = make_t()
    f = make_mountain_range_f(np.linspace(0.01, 0.99, 30))
    f.step_days = 0.03 / 2.0  # Half of the expected period
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    assert len(t) == len(y) == 30

def test_finding_maxima_near_edges():
    t0, t1 = make_t()
    f = make_mountain_range_f([bump, 1.0 - bump])
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    assert is_close(t.tt, (bump, 1.0 - bump))
    assert is_close(y, 0.0)

def test_finding_no_maxima_at_all_but_having_near_misses():
    t0, t1 = make_t()
    f = make_mountain_range_f([-bump, 1.0 + bump])
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    assert list(t.tt) == []
    assert list(y) == []

def test_finding_no_maxima_at_all_with_no_near_misses():
    t0, t1 = make_t()
    f = make_mountain_range_f([-100, 101])
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    assert list(t.tt) == []
    assert list(y) == []

def test_that_we_ignore_maxima_slightly_beyond_range():
    t0, t1 = make_t()
    f = make_mountain_range_f([-bump, 1.0 + bump])
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    assert len(t.tt) == 0
    assert len(y) == 0

def test_we_only_get_one_result_for_a_jagged_maximum():
    t0, t1 = make_t()
    almost = 0.49 * epsilon
    f = make_mountain_range_f([0.5 - almost, 0.5 + almost])
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    assert len(t.tt) == len(y) == 1

def test_we_get_two_results_for_barely_separate_maxima():
    t0, t1 = make_t()
    enough = 1.51 * epsilon
    f = make_mountain_range_f([0.5 - enough, 0.5 + enough])
    t, y = find_maxima(t0, t1, f, epsilon, 12)
    print(list(t.tt))
    assert len(t.tt) == len(y) == 2
