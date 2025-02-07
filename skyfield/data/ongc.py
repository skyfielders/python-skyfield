from math import pi

PANDAS_MESSAGE = """Skyfield needs Pandas to load the OpenNGC catalog

To load the OpenNGC star catalog, Skyfield needs the Pandas data
analysis toolkit.  Try installing it using your usual Python package
installer, like "pip install pandas" or "conda install pandas".
"""

def load_dataframe(in_dataframe):
    """Convert a dataframe from pyongc.data for import in skyfield."""
    try:
        from pandas import read_csv
    except ImportError:
        raise ImportError(PANDAS_MESSAGE)

    df = in_dataframe[['name', 'ra', 'dec', 'parallax', 'pmra', 'pmdec']].copy()

    df = df.assign(
        ra_hours = df['ra'] * 180 / pi / 15.0,
        dec_degrees = df['dec'] * 180 / pi,
        epoch_year = 2000,
    )

    df.rename(columns = {'parallax':'parallax_mas',
                         'pmra':'ra_mas_per_year',
                         'pmdec':'dec_mas_per_year'}, inplace = True)

    return df.set_index('name')
