import gzip
from skyfield.functions import to_spherical
from skyfield.starlib import Star
from skyfield.timelib import T0
from skyfield.units import Angle

days = T0 - 2448349.0625
URL = 'http://cdsarc.u-strasbg.fr/ftp/cats/I/239/hip_main.dat.gz'
url = URL  # old name, in case anyone used it

def parse(line):
    "DEPRECATED; see :func:`~skyfield.data.hipparcos.load_dataframe() instead."
    # See ftp://cdsarc.u-strasbg.fr/cats/I/239/ReadMe
    star = Star(
        ra=Angle(degrees=float(line[51:63])),
        dec=Angle(degrees=float(line[64:76])),
        ra_mas_per_year=float(line[87:95]),
        dec_mas_per_year=float(line[96:104]),
        parallax_mas=float(line[79:86]),
        names=[('HIP', int(line[8:14]))],
        )
    star._position_au += star._velocity_au_per_d * days
    distance, dec, ra = to_spherical(star._position_au)
    star.ra = Angle(radians=ra, preference='hours')
    star.dec = Angle(radians=dec)
    return star

def load(match_function):
    "DEPRECATED; see :func:`~skyfield.data.hipparcos.load_dataframe() instead."
    from skyfield import api

    with api.load.open(url) as f:
        for line in gzip.GzipFile(fileobj=f):
            if match_function(line):
                yield parse(line)

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

def load_dataframe(fobj, compression='gzip'):
    """Given an open file for `hip_main.dat.gz`, return a parsed dataframe.

    If your copy of ``hip_main.dat`` has already been unzipped, pass the
    optional argument ``compression=None``.

    """
    try:
        from pandas import read_csv
    except ImportError:
        raise ImportError(PANDAS_MESSAGE)

    df = read_csv(
        fobj, sep='|', compression=compression, names=_COLUMN_NAMES,
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

def get(which):
    "DEPRECATED; see :func:`~skyfield.data.hipparcos.load_dataframe() instead."
    if isinstance(which, str):
        pattern = ('H|      %6s' % which).encode('ascii')
        for star in load(lambda line: line.startswith(pattern)):
            return star
    else:
        patterns = set(id.encode('ascii').rjust(6) for id in which)
        return list(load(lambda line: line[8:14] in patterns))
