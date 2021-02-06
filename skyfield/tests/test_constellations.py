from skyfield.api import load_constellation_map, position_of_radec
from skyfield.timelib import Time, julian_date_of_besselian_epoch

def test_constellations():
    lookup = load_constellation_map()

    assert lookup(position_of_radec(0, 0)) == 'Psc'
    assert lookup(position_of_radec(360, 90)) == 'UMi'

    # (4.65, 0) used to be in Orion, but, precession
    assert lookup(position_of_radec(4.65, 0)) == 'Eri'
    assert lookup(position_of_radec(4.75, 0.3)) == 'Ori'

    # exploit extend() bug, issue 547
    B1875 = Time(None, julian_date_of_besselian_epoch(1875))
    assert lookup(position_of_radec(18.1747, 39., epoch=B1875)) == 'Her'
    assert lookup(position_of_radec(18.1753, 39., epoch=B1875)) == 'Lyr'

    # Test vectorization.
    assert list(
        lookup(position_of_radec([4.65, 4.75], [0, 0.3]))
    ) == ['Eri', 'Ori']
