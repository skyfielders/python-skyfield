import gzip
from skyfield.data import cache as default_cache
from skyfield.functions import to_polar
from skyfield.starlib import Star
from skyfield.timelib import T0
from skyfield.units import Angle

days = T0 - 2448349.0625
url = 'ftp://cdsarc.u-strasbg.fr/cats/I/239/hip_main.dat.gz'

def parse(line):
    # See ftp://cdsarc.u-strasbg.fr/cats/I/239/ReadMe
    star = Star(
        ra=Angle(degrees=float(line[51:63])),
        dec=Angle(degrees=float(line[64:76])),
        ra_mas_per_year=float(line[87:95]),
        dec_mas_per_year=float(line[96:104]),
        parallax_mas=float(line[79:86]),
        names=[('HIP', int(line[8:14]))],
        )
    star._position_AU += star._velocity_AU_per_d * days
    distance, dec, ra = to_polar(star._position_AU)
    star.ra = Angle(radians=ra, preference='hours')
    star.dec = Angle(radians=dec)
    return star

def load(is_match, cache=default_cache):
    with cache.open_url(url, days_old=9999) as f:
        for line in gzip.open(f):
            id = line[8:14]
            if is_match(id):
                yield parse(line)

def get(which, cache=default_cache):
    #print(repr(which).rjust(6))
    if isinstance(which, str):
        is_match = which.rjust(6).encode('ascii').__eq__
        for star in load(is_match):
            return star
    else:
        id_set = set(str(id).encode('ascii').rjust(6) for id in which)
        is_match = id_set.__contains__
        return list(load(is_match))
