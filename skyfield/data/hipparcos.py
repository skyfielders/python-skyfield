# This URL worked until September 2020:
#
# URL = 'http://cdsarc.u-strasbg.fr/ftp/cats/I/239/hip_main.dat.gz'
#
# Then someone at VizieR apparently ran `gunzip` on the file, breaking
# the existing URL.  The fastest fix is for us to switch to:
#
# URL = 'https://cdsarc.u-strasbg.fr/ftp/cats/I/239/hip_main.dat'
#
# Then in September 2022 the site's certificate broke, and an engineer
# at unistra.fr confirmed that 'the domain "astro.unistra.fr" as well as
# the domain "u-strasbg.fr" are obsoleted by "cds.unistra.fr".'  Thus:

URL = 'https://cdsarc.cds.unistra.fr/ftp/cats/I/239/hip_main.dat'

# But what if someone runs `gzip` on the file again?  Then the new URL
# will break like the old one did.  It appears that VizieR makes no
# guarantee that raw catalog URLs are stable, and that we need to switch
# to one of their catalog generation URLs, which produce text in a new
# format that Skyfield will have to learn.  Discussion at:
# https://github.com/skyfielders/python-skyfield/issues/454

url = URL  # old name, in case anyone used it

PANDAS_MESSAGE = """Skyfield needs Pandas to load the Hipparcos catalog

To load the Hipparcos star catalog, Skyfield needs the Pandas data
analysis toolkit.  Try installing it using your usual Python package
installer, like "pip install pandas" or "conda install pandas".
"""

_COLUMN_NAMES = (
    'Catalog', 'HIP', 'Proxy', 'RAhms', 'DEdms', 'Vmag',
    'VarFlag', 'r_Vmag', 'RAdeg', 'DEdeg', 'AstroRef', 'Plx', 'pmRA',
    'pmDE', 'e_RAdeg', 'e_DEdeg', 'e_Plx', 'e_pmRA', 'e_pmDE', 'DE:RA',
    'Plx:RA', 'Plx:DE', 'pmRA:RA', 'pmRA:DE', 'pmRA:Plx', 'pmDE:RA',
    'pmDE:DE', 'pmDE:Plx', 'pmDE:pmRA', 'F1', 'F2', '---', 'BTmag',
    'e_BTmag', 'VTmag', 'e_VTmag', 'm_BTmag', 'B-V', 'e_B-V', 'r_B-V',
    'V-I', 'e_V-I', 'r_V-I', 'CombMag', 'Hpmag', 'e_Hpmag', 'Hpscat',
    'o_Hpmag', 'm_Hpmag', 'Hpmax', 'HPmin', 'Period', 'HvarType',
    'moreVar', 'morePhoto', 'CCDM', 'n_CCDM', 'Nsys', 'Ncomp',
    'MultFlag', 'Source', 'Qual', 'm_HIP', 'theta', 'rho', 'e_rho',
    'dHp', 'e_dHp', 'Survey', 'Chart', 'Notes', 'HD', 'BD', 'CoD',
    'CPD', '(V-I)red', 'SpType', 'r_SpType',
)

def load_dataframe(fobj):
    """Given an open file for ``hip_main.dat``, return a parsed dataframe.

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
        usecols=['HIP', 'Vmag', 'RAdeg', 'DEdeg', 'Plx', 'pmRA', 'pmDE'],
        na_values=['     ', '       ', '        ', '            '],
    )
    df.columns = (
        'hip', 'magnitude', 'ra_degrees', 'dec_degrees',
        'parallax_mas', 'ra_mas_per_year', 'dec_mas_per_year',
    )
    df = df.assign(
        ra_hours = df['ra_degrees'] / 15.0,
        epoch_year = 1991.25,
    )
    return df.set_index('hip')
