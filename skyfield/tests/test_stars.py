from skyfield import api
from skyfield.data.hipparcos import load_dataframe

def test_dataframe():
    with api.load.open('hip_main.dat.gz') as f:
        df = load_dataframe(f)
    star = api.Star.from_dataframe(df)
    assert repr(star) == 'Star(ra shape=9933, dec shape=9933, ra_mas_per_year shape=9933, dec_mas_per_year shape=9933, parallax_mas shape=9933, epoch shape=9933)'
