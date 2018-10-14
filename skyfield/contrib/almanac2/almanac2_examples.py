from skyfield.api import Loader, Topos, EarthSatellite, Star
from almanac2 import seasons, moon_phases, meridian_transits, culminations, twilights, risings_settings

load = Loader(r'C:\Users\Josh\Scripts\Skyfield_Data')
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


# Equinoxes and Solstices
season_times, lons = seasons(earth, ts.utc(2000), ts.utc(2026))

# They can be separated like this:
march_equinoxes = season_times[lons.degrees==0]
june_solstices = season_times[lons.degrees==90]
sept_equinoxs = season_times[lons.degrees==180]
dec_solstices = season_times[lons.degrees==270]


# Moon Phases
phase_times, lon_diffs = moon_phases(moon, ts.utc(2018, 1), ts.utc(2018, 2))

# They can be separated like this:
new_moons = phase_times[lon_diffs.degrees==0]
first_quarters = phase_times[lon_diffs.degrees==90]
full_moons = phase_times[lon_diffs.degrees==180]
last_quarters = phase_times[lon_diffs.degrees==270]


# Rise and Set Times
for name, body in zip(['sun', 'mars', "Barnard's Star"], [sun, mars, barnard]):
    times, kinds = risings_settings(greenwich, body, ts.utc(2018, 1, 1), ts.utc(2018, 1, 2))
    print (name, ' rises at:', times[kinds=='rise'].utc_jpl)
