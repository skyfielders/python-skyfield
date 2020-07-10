from skyfield import api
from skyfield.data.hipparcos import load_dataframe

def test_dataframe():
    with api.load.open('hip_main.dat.gz') as f:
        df = load_dataframe(f)
    star = api.Star.from_dataframe(df.iloc[:214])
    assert repr(star) == 'Star(ra shape=214, dec shape=214, ra_mas_per_year shape=214, dec_mas_per_year shape=214, parallax_mas shape=214, epoch shape=214)'
