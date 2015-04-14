import gzip
from skyfield.data import cache as default_cache
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

def load(match_function, cache=default_cache):
    """Yield the Hipparcos stars for which `match_function(line)` is true."""
    with cache.open_url(url, days_old=9999) as f:
        for line in gzip.GzipFile(fileobj=f):
            if match_function(line):
                yield parse(line)

def get(which, cache=default_cache):
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
