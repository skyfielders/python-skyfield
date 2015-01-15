"""Tests for routines from the functions module."""

from skyfield.constants import tau
from skyfield.functions import from_polar

def test_right():
    dx, dy, dz = from_polar(1, 0, 0) - [1, 0, 0]
    assert abs(dx) < 1e-15
    assert abs(dy) < 1e-15
    assert abs(dz) < 1e-15

def test_down():
    dx, dy, dz = from_polar(1, 0, 0.75 * tau) - [0, -1, 0]
    assert abs(dx) < 1e-15
    assert abs(dy) < 1e-15
    assert abs(dz) < 1e-15

def test_up():
    dx, dy, dz = from_polar(1, 0.25 * tau, 0) - [0, 0, 1]
    assert abs(dx) < 1e-15
    assert abs(dy) < 1e-15
    assert abs(dz) < 1e-15

def test_left_up():
    sqrt2 = 0.5 ** 0.5
    dx, dy, dz = from_polar(1, 0.125 * tau, 0.5 * tau) - [-sqrt2, 0, sqrt2]
    assert abs(dx) < 1e-15
    assert abs(dy) < 1e-15
    assert abs(dz) < 1e-15
