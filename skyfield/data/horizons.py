"""Physical data for the planets from the JPL HORIZONS system.

To rebuild this data, consult the following IPython Notebook:

https://github.com/brandon-rhodes/astronomy-notebooks/blob/master/Utils-HORIZONS-data.ipynb

"""
from skyfield.units import Distance

radii_km = [
    ('Sun', 695500.0),
    ('Mercury', 2440.0),
    ('Venus', 6051.8),
    ('Earth', 6371.01),
    ('Mars', 3389.9),
    ('Jupiter', 69911.0),
    ('Saturn', 58232.0),
    ('Uranus', 25362.0),
    ('Neptune', 24624.0),
    ('134340 Pluto', 1195.0),
    ]

def festoon_ephemeris(ephemeris):
    for name, radius_km in radii_km:
        name = name.lower().split()[-1]
        getattr(ephemeris, name).radius = Distance(km=radius_km)
