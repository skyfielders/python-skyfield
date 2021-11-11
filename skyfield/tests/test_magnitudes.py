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
    assert ['%.3f' % m for m in magnitudes] == [
        '2.393', '-4.278', '1.592', '-2.693', '0.508', '5.701', '7.690',
    ]
