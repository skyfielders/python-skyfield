from skyfield.api import (
    load_constellation_map, load_constellation_names, position_of_radec,
)
from skyfield.functions import load_bundled_npy
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

    # Verify we are not using the paper's original table, which had an error.
    star = position_of_radec(16.323165, -18.848203, epoch=B1875)
    assert lookup(star) == 'Sco'

    # Test vectorization.
    assert list(
        lookup(position_of_radec([4.65, 4.75], [0, 0.3]))
    ) == ['Eri', 'Ori']

def test_that_constellation_abbreviations_match():
    # This test should ideally not have to load the underlying data
    # source that is behind load_constellation_map(), but for the moment
    # I do not see a quick alternative.
    arrays = load_bundled_npy('constellations.npz')
    abbrevs1 = set(arrays['indexed_abbreviations'])

    abbrevs2 = {abbrev for abbrev, name in load_constellation_names()}
    assert abbrevs1 == abbrevs2
