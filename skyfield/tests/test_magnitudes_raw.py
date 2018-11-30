
from skyfield import magnitudelib as m

def test_mercury_magnitude_function():
    assert abs(-2.477 - m.mercury_magnitude(0.310295423552, 1.32182643625754, 1.1677)) < 0.0005
    assert abs(0.181 - m.mercury_magnitude(0.413629222334, 0.92644808718613, 90.1662)) < 0.0005
    assert abs(7.167 - m.mercury_magnitude(0.448947624811, 0.56004973217883, 178.7284)) < 0.0005

def test_venus_magnitude_function():
    assert abs(-3.917 - m.venus_magnitude(0.722722540169, 1.71607489554051, 1.3232)) < 0.0005
    assert abs(-4.916 - m.venus_magnitude(0.721480714554, 0.37762511206278, 124.1348)) < 0.0005
    assert abs(-3.090 - m.venus_magnitude(0.726166592736, 0.28889582420642, 179.1845)) < 0.0005
