import numpy as np
from skyfield.api import load
from skyfield.constants import DAY_S
from skyfield.elementslib import osculating_elements_of
from skyfield.positionlib import ICRF
from skyfield.units import Distance, Velocity


def _hyperbolic_position(H):
    # Parameterize a hyperbola directly, independently of ele_to_vec().
    a, e, mu = 1e6, 1.3, 398600.4418
    n = np.sqrt(mu / a**3)
    denominator = e * np.cosh(H) - 1.0
    r = a * np.array([e - np.cosh(H),
                      np.sqrt(e**2 - 1.0) * np.sinh(H),
                      np.zeros_like(H)])
    v = a * n / denominator * np.array([
        -np.sinh(H), np.sqrt(e**2 - 1.0) * np.cosh(H), np.zeros_like(H)])
    ts = load.timescale()
    t = ts.tdb_jd(2451545.0, (e * np.sinh(H) - H) / n / DAY_S)
    return ICRF(Distance(km=r).au, Velocity(km_per_s=v).au_per_d, t), mu


def test_hyperbolic_anomaly_scalar():
    for H in [-4.0, -np.pi, -1.0, 0.0, 1.0, np.pi, 4.0]:
        position, mu = _hyperbolic_position(H)
        elements = osculating_elements_of(position, gm_km3_s2=mu)
        assert abs(elements.eccentric_anomaly.radians - H) < 1e-12
        assert abs(elements.periapsis_time.tdb - 2451545.0) < 1e-8


def test_hyperbolic_anomaly_array():
    H = np.array([-4.0, -np.pi, -1.0, 0.0, 1.0, np.pi, 4.0])
    position, mu = _hyperbolic_position(H)
    elements = osculating_elements_of(position, gm_km3_s2=mu)
    assert np.max(np.abs(elements.eccentric_anomaly.radians - H)) < 1e-12
    assert np.max(np.abs(elements.periapsis_time.tdb - 2451545.0)) < 1e-8
