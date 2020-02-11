from numpy import pi, seterr, linspace
from skyfield.keplerlib import propagate
from skyfield.elementslib import OsculatingElements
from skyfield.units import Velocity, Distance
from skyfield.tests.test_elementslib import compare, ele_to_vec
from skyfield.api import load
import os

seterr(all='raise')

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

def check_orbit(p, e, i, Om, w, v, 
                p_eps=None, e_eps=None, i_eps=None, Om_eps=None, w_eps=None, v_eps=None):
    pos0, vel0 = ele_to_vec(p, e, i, Om, w, v, mu)
    
    pos1, vel1 = propagate(pos0, vel0, 0, times, mu)
    ele = OsculatingElements(Distance(km=pos1), Velocity(km_per_s=vel1), dummy_time, mu)
    
    if p_eps: compare(p, ele.semi_latus_rectum.km, p_eps)
    if e_eps: compare(e, ele.eccentricity, e_eps)
    if i_eps: compare(i, ele.inclination.radians, i_eps, mod=True)
    if Om_eps: compare(Om, ele.longitude_of_ascending_node.radians, Om_eps, mod=True)
    if w_eps: compare(w, ele.argument_of_periapsis.radians, w_eps, mod=True)
    if v_eps: compare(v, ele.true_anomaly.radians, v_eps, mod=True)
    

times = linspace(-1e11, 1e11, 1001) # -3170 years to +3170 years, including 0
mu = 403503.2355022598
dummy_time = load.timescale().utc(2018)

def test_circular():
    check_orbit(p=300000, e=0, i=.5, Om=1, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)


def test_circular_equatorial():
    check_orbit(p=300000, e=0, i=0, Om=0, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15)


def test_circular_polar():
    check_orbit(p=300000, e=0, i=pi/2, Om=1, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)
        
        
def test_elliptical():
    check_orbit(p=300000, e=.3, i=1, Om=0, w=4, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)
        

def test_elliptical_equatorial():
    check_orbit(p=300000, e=.3, i=0, Om=0, w=1, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)
        

def test_elliptical_polar():
    check_orbit(p=300000, e=.2, i=pi/2, Om=1, w=2, v=3,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-8)
        

def test_parabolic():
    check_orbit(p=300000, e=1, i=1, Om=0, w=4, v=3,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-13, Om_eps=1e-13, w_eps=1e-13)
        

def test_parabolic_equatorial():
    check_orbit(p=300000, e=1, i=0, Om=0, w=1, v=2,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-13)
        

def test_parabolic_polar():
    check_orbit(p=300000, e=1, i=pi/2, Om=1, w=2, v=3,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-14, Om_eps=1e-13, w_eps=1e-13)


def test_hyperbolic():
    check_orbit(p=300000, e=1.3, i=1, Om=0, w=4, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-10, Om_eps=1e-10, w_eps=1e-6)        


def test_hyperbolic_equatorial():
    check_orbit(p=300000, e=1.3, i=0, Om=0, w=1, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-6)


def test_hyperbolic_polar():
    check_orbit(p=300000, e=1.3, i=pi/2, Om=1, w=2, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-10, Om_eps=1e-10, w_eps=1e-6)
