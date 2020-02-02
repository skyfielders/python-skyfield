"""Low-level tests of the almanac search routines."""

from numpy import where, sin

from skyfield.api import load
from skyfield.constants import tau
from skyfield.searchlib import find_discrete, _find_maxima as find_maxima

bump = 1e-5
epsilon = 1e-10

def make_t():
    ts = load.timescale(builtin=True)
    t0 = ts.tt_jd(0)
    t1 = ts.tt_jd(1)
    return t0, t1

def make_f(offset):
    """Make a sine wave from 0 to 1, with `t` offset by `offset` days."""
    def f(t):
        return sin((t.tt + offset) * tau) >= 0.0
    f.rough_period = 1.0
    return f

def is_close(value, expected):
    return (abs(value - expected) < epsilon).all()

def test_find_discrete_near_left_edge():
    t0, t1 = make_t()
    f = make_f(-bump)  # cross zero barely past t0
    t, y = find_discrete(t0, t1, f, epsilon)
    assert is_close(t.tt, (bump, 0.5 + bump))
    assert list(y) == [1, 0]

def test_find_discrete_near_right_edge():
    t0, t1 = make_t()
    f = make_f(bump)  # cross zero almost at the end of the range
    t, y = find_discrete(t0, t1, f, epsilon)
    assert is_close(t.tt, (0.5 - bump, 1.0 - bump))
    assert list(y) == [0, 1]

def test_find_discrete_with_a_barely_detectable_jag_right_at_zero():
    t0, t1 = make_t()
    def f(t):
        n = t.tt
        n = where(n < 0.5, n + 3.1 * epsilon, n - 3.1 * epsilon)
        return sin(n * tau) >= 0.0
    f.rough_period = 1.0
    t, y = find_discrete(t0, t1, f, epsilon)
    assert is_close(t.tt, (0.5 - 3.1 * epsilon, 0.5, 0.5 + 3.1 * epsilon))
    assert list(y) == [0, 1, 0]

def test_find_discrete_with_a_sub_epsilon_jag_right_at_zero():
    t0, t1 = make_t()
    def f(t):
        n = t.tt
        n = where(n < 0.5, n + epsilon * 0.99, n - epsilon * 0.99)
        return sin(n * tau) >= 0.0
    f.rough_period = 1.0

    # We hard-code num=12, just in case the default ever changes to
    # another value that might not trigger the symptom.
    t, y = find_discrete(t0, t1, f, epsilon, 12)

    assert is_close(t.tt, (0.5,))
    assert list(y) == [0]
