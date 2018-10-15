from skyfield.api import Loader, Topos, EarthSatellite, Star
from almanac2 import seasons, moon_phases, meridian_transits, culminations, twilights, risings_settings

# Put your data directory here before running this file:
load = Loader('')

ts = load.timescale()

ephem = load('de430.bsp')
earth = ephem['earth']
sun = ephem['sun']
moon = ephem['moon']
mars = ephem['mars barycenter']

greenwich = earth + Topos(latitude_degrees=(51, 28, 40), longitude_degrees=(0, 0, -5))

iss_tle = """\
1 25544U 98067A   18161.85073725  .00003008  00000-0  52601-4 0  9993
2 25544  51.6418  50.3007 0003338 171.6979 280.7366 15.54163173117534
"""
ISS = EarthSatellite(*iss_tle.splitlines())

sirius = Star(ra_hours=(6, 45, 8.91728), dec_degrees=(-16, 42, 58.0171))

t0 = ts.utc(2017, 1, 1)
t1 = ts.utc(2017, 1, 8)


# Equinoxes and Solstices
season_times, lons = seasons(earth, ts.utc(2000), ts.utc(2026))

# Moon Phases
phase_times, lon_diffs = moon_phases(moon, ts.utc(2018, 1), ts.utc(2018, 2))

# Rise and Set Times
rise_set_times, rise_or_set = risings_settings(greenwich, sun, t0, t1)

# Culminations
culmination_times, kinds = culminations(greenwich, sun, t0, t1)

# Meridian Transits
transit_times, hour_angles = meridian_transits(greenwich, mars, t0, t1)

# Twilights
naut_twilight_times, am_pm = twilights(greenwich, sun, t0, t1, kind='nautical')

# You can submit different object types to the same function:
for name, body in zip(['Sun', 'Mars', 'Sirius', 'ISS'], [sun, mars, sirius, ISS]):
    times, kinds = risings_settings(greenwich, body, t0, t1)
    print (name, 'first rises at:', times[kinds=='rise'][0].utc_jpl())
