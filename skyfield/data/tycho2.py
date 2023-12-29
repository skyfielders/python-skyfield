URL = 'https://cdsarc.u-strasbg.fr/ftp/cats/I/239/tyc_main.dat'

PANDAS_MESSAGE = """Skyfield needs Pandas to load the Tycho2 catalog

To load the Tycho2 star catalog, Skyfield needs the Pandas data
analysis toolkit.  Try installing it using your usual Python package
installer, like "pip install pandas" or "conda install pandas".
"""

_COLUMN_NAMES = (
    'Catalog', 'TYC', 'Proxy', 'RAhms', 'DEdms', 'Vmag', '---',
    'r_Vmag', 'RAdeg', 'DEdeg', 'AstroRef', 'Plx', 'pmRA',
    'pmDE', 'e_RAdeg', 'e_DEdeg', 'e_Plx', 'e_pmRA', 'e_pmDE', 'DE:RA',
    'Plx:RA', 'Plx:DE', 'pmRA:RA', 'pmRA:DE', 'pmRA:Plx', 'pmDE:RA',
    'pmDE:DE', 'pmDE:Plx', 'pmDE:pmRA', 'Nastro', 'F2', 'HIP', 'BTmag',
    'e_BTmag', 'VTmag', 'e_VTmag', 'r_BTmag', 'B-V', 'e_B-V', '---(2)',
    'Q', 'Fs', 'Source', 'Nphoto', 'VTscat', 'VTmax', 'VTmin', 'Var',
    'VarFlag', 'MultFlag', 'morePhoto', 'm_HIP', 'PPM', 'HD', 'BD', 'CoD',
    'CPD', 'Remark',
)

def load_dataframe(fobj):
    """Given an open file for ``tyc_main.dat``, return a parsed dataframe.

    If the file is gzipped, it will be automatically uncompressed.

    """
    try:
        from pandas import read_csv
    except ImportError:
        raise ImportError(PANDAS_MESSAGE)

    fobj.seek(0)
    magic = fobj.read(2)
    compression = 'gzip' if (magic == b'\x1f\x8b') else None
    fobj.seek(0)

    df = read_csv(
        fobj, sep='|', names=_COLUMN_NAMES, compression=compression,
        usecols=['TYC', 'Vmag', 'RAdeg', 'DEdeg', 'Plx', 'pmRA', 'pmDE'],
        na_values=['     ', '       ', '        '],
    )
    df.columns = (
        'tyc', 'magnitude', 'ra_degrees', 'dec_degrees',
        'parallax_mas', 'ra_mas_per_year', 'dec_mas_per_year',
    )
    df = df.assign(
        ra_hours = df['ra_degrees'] / 15.0,
        epoch_year = 1991.25,
    )
    return df.set_index('tyc')
