from skyfield import api
from skyfield.data.tycho2 import load_dataframe

def test_dataframe():
    with api.load.open('tyc_main_head.dat') as f:
        df = load_dataframe(f)
    star = api.Star.from_dataframe(df)
    assert repr(star) == 'Star(ra shape=10, dec shape=10, ra_mas_per_year shape=10, dec_mas_per_year shape=10, parallax_mas shape=10, epoch shape=10)'
