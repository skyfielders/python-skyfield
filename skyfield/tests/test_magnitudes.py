from skyfield.api import load
from skyfield.magnitudelib import planetary_magnitude

def test_magnitudes():
    ts = load.timescale()
    #t = ts.utc(1995, 5, 22)  # Rings edge-on from Earth.
    t = ts.utc(2021, 10, 4)
    eph = load('de421.bsp')
    names = [
        'mercury', 'venus', 'mars', 'jupiter barycenter',
        'saturn barycenter', 'uranus barycenter', 'neptune barycenter',
    ]
    e = eph['earth'].at(t)
    positions = [e.observe(eph[name]) for name in names]
    magnitudes = [planetary_magnitude(position) for position in positions]
    assert [f'{m:.1f}' for m in magnitudes] == [
        '2.4', '-4.3', '1.6', '-2.7', '0.5', '5.7', '7.7',
    ]
