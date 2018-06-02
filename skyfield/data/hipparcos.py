import gzip
from skyfield.functions import to_polar
from skyfield.starlib import Star
from skyfield.timelib import T0
from skyfield.units import Angle

days = T0 - 2448349.0625
url = 'ftp://cdsarc.u-strasbg.fr/cats/I/239/hip_main.dat.gz'

def parse(line):
    """Return a `Star` build by parsing a Hipparcos catalog entry `line`."""
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
    distance, dec, ra = to_polar(star._position_au)
    star.ra = Angle(radians=ra, preference='hours')
    star.dec = Angle(radians=dec)
    return star

def load(match_function):
    """Yield the Hipparcos stars for which `match_function(line)` is true."""
    from skyfield import api

    with api.load.open(url) as f:
        for line in gzip.GzipFile(fileobj=f):
            if match_function(line):
                yield parse(line)

def load_dataframe(path):
    from pandas import read_fwf
    # from csv import reader

    # with api.load.open(url) as f:
    #     df = read_csv(f, sep='|', compression='gzip')

    from time import time
    t0 = time()

    # with api.load.open(url) as f:
    #     #for line in reader(gzip.open(f, 'rt', encoding='ascii')):
    #     for line in gzip.open(f, 'rt', encoding='ascii'):
    #         # assert len(line) >= 43, repr(line)
    #         assert line[43] == '.', (repr(line), repr(line[43]))
    #         if match_function(line):
    #             yield parse(line)

    names, colspecs = zip(
        ('hip', (2, 14)),
        ('magnitude', (41, 46)),
        ('ra_degrees', (51, 63)),
        ('dec_degrees', (64, 76)),
    )

    #with api.load.open(url) as f:
    with open(path, 'rb') as f:
        with gzip.open(f, 'rt', encoding='ascii') as g:
        # with gzip.open(f, 'rb') as g:
            #print(g)
            df = read_fwf(g, names=names, colspecs=colspecs)

    #'H|         104| |00 01 18.42|-31 53 58.6| 9.'
    #print(df.head(1))
    print(time() - t0)
    df = df.set_index('hip')
    df['ra_hours'] = df['ra_degrees'] / 15.0
    # print(df.head())
    # print(df.info())
    # print(df.magnitude.isnull().sum())
    # print(df.ra_hours.isnull().sum())
    # print(df.dec_degrees.isnull().sum())
        # for line in gzip.GzipFile(fileobj=f):
        #     if match_function(line):
        #         yield parse(line)
    return df

def get(which):
    """Return a single star, or a list of stars, from the Hipparcos catalog.

    A call like `get('54061')` returns a single `Star` object, while
    `get(['54061', '53910'])` returns a list of stars.

    """
    if isinstance(which, str):
        pattern = ('H|      %6s' % which).encode('ascii')
        for star in load(lambda line: line.startswith(pattern)):
            return star
    else:
        patterns = set(id.encode('ascii').rjust(6) for id in which)
        return list(load(lambda line: line[8:14] in patterns))
